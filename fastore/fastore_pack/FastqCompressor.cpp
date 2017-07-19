/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM and QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../fastore_bin/Globals.h"

#include <string.h>
#include <queue>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "FastqCompressor.h"
#include "CompressedBlockData.h"

#include "../fastore_bin/FastqPacker.h"
#include "../rle/rle.h"
#include "../ppmd/PPMd.h"
#include "quantizer.h"


#include <stdio.h>


// IFastqCompressor
//
IStoreBase::IStoreBase(const CompressorParams &params_, const CompressorAuxParams &auxParams_)
	:	params(params_)
	,	auxParams(auxParams_)
	,	dryFastqBuffer(NULL)
	,	dryFastqWriter(NULL)
	,	currentRecordIdx(0)
{
	if (auxParams_.dry_run)
	{
		dryFastqBuffer = new Buffer(FastqChunk::DefaultBufferSize);
		dryFastqWriter = new BitMemoryWriter(*dryFastqBuffer);
	}
}

IStoreBase::~IStoreBase()
{
	if (dryFastqWriter != NULL)
		delete dryFastqWriter;

	if (dryFastqBuffer != NULL)
		delete dryFastqBuffer;
}


void IStoreBase::StoreRawHeader(const BaseBlockHeader& header_, BitMemoryWriter& writer_)
{
	writer_.Put4Bytes(header_.minimizerId);
	writer_.Put8Bytes(header_.recordsCount);			// this goes to META
	writer_.PutByte(header_.recMinLen);					// -- * --
	writer_.PutByte(header_.recMaxLen);					// -- * --
	writer_.Put8Bytes(header_.rawDnaStreamSize);

	writer_.Put8Bytes(header_.footerOffset);
	writer_.Put4Bytes(header_.footerSize);

	if (params.archType.readsHaveHeaders)
		writer_.Put8Bytes(header_.rawIdStreamSize);
}

void IStoreBase::ReadRawHeader(BaseBlockHeader& header_, BitMemoryReader& reader_)
{
	header_.minimizerId = reader_.Get4Bytes();
	header_.recordsCount = reader_.Get8Bytes();
	header_.recMinLen = reader_.GetByte();
	header_.recMaxLen = reader_.GetByte();
	header_.rawDnaStreamSize = reader_.Get8Bytes();

	header_.footerOffset = reader_.Get8Bytes();
	header_.footerSize = reader_.Get4Bytes();

	if (params.archType.readsHaveHeaders)
	{
		header_.rawIdStreamSize = reader_.Get8Bytes();
		ASSERT(header_.rawIdStreamSize > 0);
		ASSERT(reader_.Position() == BaseBlockHeader::Size);
	}
	else
	{
		header_.rawIdStreamSize = 0;
		ASSERT(reader_.Position() == BaseBlockHeader::Size - sizeof(header_.rawIdStreamSize));
	}

	ASSERT(header_.footerOffset + header_.footerSize <= reader_.Size());
	ASSERT(header_.recordsCount > 0);
	ASSERT(header_.recMaxLen > 0);
	ASSERT(header_.recMaxLen >= header_.recMinLen);
	ASSERT(header_.rawDnaStreamSize > 0);
}


void IStoreBase::StoreRawFooter(const BaseBlockFooter& footer_, BitMemoryWriter& writer_)
{
	writer_.PutByte(footer_.sampleValue);
}


void IStoreBase::ReadRawFooter(BaseBlockFooter& footer_, BitMemoryReader& reader_)
{
	footer_.sampleValue = reader_.GetByte();
}


void IStoreBase::CompressBuffer(PpmdEncoder& encoder_, const DataChunk& inChunk_, uint64 inSize_,
								   DataChunk& outChunk_, uint64& outSize_, uint64 outOffset_)
{
	ASSERT(inSize_ > 0);
	outSize_ = 0;

	if (outOffset_ + inSize_ > outChunk_.data.Size())
	{
		outChunk_.data.Extend(outOffset_ + inSize_, true);
	}

	byte* inMem = inChunk_.data.Pointer();
	byte* outMem =  outChunk_.data.Pointer() + outOffset_;
	uint64_t outSize = outChunk_.data.Size() - outOffset_;

	bool r = encoder_.EncodeNextMember(inMem, inSize_, outMem, outSize);
	ASSERT(r);
	ASSERT(outSize > 0);

	outSize_ = outSize;
}

void IStoreBase::DecompressBuffer(PpmdDecoder& decoder_, DataChunk& outChunk_, uint64& outSize_,
									 const DataChunk& inChunk_, uint64 inSize_, uint64 inOffset_)
{
	ASSERT(inSize_ > 0);

	byte* inMem = inChunk_.data.Pointer() + inOffset_;

	byte* outMem = outChunk_.data.Pointer();
	uint64_t outSize = outChunk_.data.Size();

	ASSERT(outSize > inSize_);

	bool r = decoder_.DecodeNextMember(inMem, (uint64_t)inSize_, outMem, outSize);
	ASSERT(r);
	ASSERT(outSize > 0);

	outSize_ = outSize;
}




IDnaStoreBase::IDnaStoreBase(const MinimizerParameters &minimizer_)
{
	// dna translation tables - those should be used in separate modules
	// to avoid code duplication
	//
	std::fill(dnaToIdx.begin(), dnaToIdx.end(), -1);
	for (uint32 i = 0; i < 5; ++i)
	{
		int32 c = minimizer_.dnaSymbolOrder[i];
		ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
		dnaToIdx[c] = i;
	}

	std::fill(idxToDna.begin(), idxToDna.end(), -1);
	for (uint32 i = 0; i < 5; ++i)
	{
		idxToDna[i] = minimizer_.dnaSymbolOrder[i];
	}
}


IQualityStoreBase::IQualityStoreBase(const QualityCompressionData& globalQuaData_)
	:	globalQuaData(globalQuaData_)
{
	// qualities translation tables -- those should be used in separate modules
	// to avoid code duplication
	//
	std::array<uint8, 64> qualityTranslationTable = {{
							0, 0, 6, 6, 6,	6, 6, 6, 6, 6,	// 0x
							15,15,15,15,15,	15,15,15,15,15,	// 1x
							22,22,22,22,22,	27,27,27,27,27,	// 2x
							33,33,33,33,33,	37,37,37,37,37,	// 3x
							40,40,40,40,40,	40,40,40,40,40,	// 4x
							40,40,40,40,40,	40,40,40,40,40,	// 5x
							40,40,40,40	}};					// 6x - 64


	idxToQua_8bin = {{0, 6, 15, 22, 27, 33, 37, 40}};

	uint8 prev = 0;
	uint32 sym = 0;

	for (uint32 i = 0; i < qualityTranslationTable.size(); ++i)
	{
		if (qualityTranslationTable[i] != prev)
		{
			prev = qualityTranslationTable[i];
			sym++;
		}
		quaToIdx_8bin[i] = sym;
	}

	ResetWellRng();
}

void IQualityStoreBase::ResetWellRng()
{
	// We need to initialize the WELL RND
	memset(&local_well, 0, sizeof(struct well_state_t));
	std::copy((byte*)&globalQuaData.well.state, (byte*)&globalQuaData.well.state + sizeof(local_well.state), (byte*)&local_well.state);
}


void IQualityStoreBase::CompressReadQuality(const FastqRecord &rec_,
											QualityEncoders &enc_,
											const CompressorParams& params_,
											const bool dryRun_,
											std::vector<byte>* dryBuffer_)
{
	ASSERT(!dryRun_ || dryBuffer_ != NULL);

	switch (params_.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			// special index to access q-scores from the end in case of rev-compl
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;
			uint32 q = rec_.qua[ii] - params_.archType.qualityOffset;
			//char c = rec_.seq[ii];

			if (dryRun_)
			{
				dryBuffer_->push_back(rec_.qua[ii]);
			}
			else
			{
				//if (c != 'N')
				//{
					//ASSERT(q != 0);
					enc_.rawCoder->PutByte(q);
				//}
			}
		}
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			// special index to access q-scores from the end in case of rev-compl
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;

			uint32 q = rec_.qua[ii];
			ASSERT(q >= params_.archType.qualityOffset && q < params_.archType.qualityOffset + 64U);
			q = q - params_.archType.qualityOffset;
			q = q >= params_.quality.binaryThreshold;

			char c = rec_.seq[ii];

			if (dryRun_)
			{
				dryBuffer_->push_back(params_.archType.qualityOffset + (q ? QualityCompressionParams::Default::MaxThresholdValue
																 : QualityCompressionParams::Default::MinThresholdValue));
			}
			else
			{
				if (c != 'N')
				{
					//ASSERT(q != 0);
					uint32 posCtx = (i * 2) / rec_.seqLen;
					enc_.binaryCoder->coder.EncodeSymbol(enc_.binaryCoder->rc, q, posCtx);
				}
			}
		}
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			// special index to access q-scores from the end in case of rev-compl
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;

			char c = rec_.seq[ii];

			uint32 q = rec_.qua[ii];
			ASSERT(q >= params_.archType.qualityOffset && q < params_.archType.qualityOffset + 64U);
			q = quaToIdx_8bin[q - params_.archType.qualityOffset];

			if (dryRun_)
			{
				dryBuffer_->push_back(params_.archType.qualityOffset + idxToQua_8bin[q]);
			}
			else
			{
				if (c != 'N')
				{
					//ASSERT(q != 0);
					uint32 posCtx = (i * 8) / rec_.seqLen;
					enc_.illu8Coder->coder.EncodeSymbol(enc_.illu8Coder->rc, q, posCtx);
				}
			}
		}
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		struct cond_quantizer_list_t *qlist = globalQuaData.codebook.qlist;
		struct quantizer_t *q;
		uint32_t idx = 0;
		uint32_t qv_hat_prev = 0;
		uint32_t qv = 0, qv_hat = 0, qv_state = 0;

		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			// special index to access q-scores from the end in case of rev-compl
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;

			qv = rec_.qua[ii] - params_.archType.qualityOffset;
			ASSERT(qv < 64);


			// choose the quantizer given the current pos and the previous quantized value
			q = choose_quantizer(qlist, &local_well, i, qv_hat_prev, &idx);

			// quantize the value
			qv_hat = q->q[qv];

			// get the state [we prefer to compress the state of the quantizer than the actual value of qv_hat]
			qv_state = get_symbol_index(q->output_alphabet, qv_hat);



			// encode the state
			if (dryRun_)
			{
				dryBuffer_->push_back(params_.archType.qualityOffset + qv_hat);
			}
			else
			{
				//enc_.qvzCoder->coder.EncodeSymbol(enc_.qvzCoder->rc, qv_state, idx);
				//ASSERT(qv != 0);

				enc_.qvzCoder->EncodeNext(qv_state, i, idx);
			}

			// store the current qv_hat for next iteration
			qv_hat_prev = qv_hat;
		}

		break;
	}
	}
}


void IQualityStoreBase::DecompressReadQuality(FastqRecord &rec_,
											  QualityDecoders &dec_,
											  const CompressorParams& params_)
{
	switch (params_.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;

			uint32 q = 0;
			//char c = rec_.seq[ii];
			//if (c != 'N')
			//{
				q = dec_.rawCoder->GetByte();
				ASSERT(q < 64);
			//}

			rec_.qua[ii] = (char)(q + params_.archType.qualityOffset);
		}
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;

			uint q = 0;
			char c = rec_.seq[ii];
			if (c != 'N')
			{
				uint32 posCtx = (i * 2) / rec_.seqLen;
				q = dec_.binaryCoder->coder.DecodeSymbol(dec_.binaryCoder->rc, posCtx);
				ASSERT(q <= 1);
			}

			rec_.qua[ii] = params_.archType.qualityOffset + (q ? QualityCompressionParams::Default::MaxThresholdValue
															 : QualityCompressionParams::Default::MinThresholdValue);
		}
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;

			uint32 q = 0;
			char c = rec_.seq[ii];

			if (c != 'N')
			{
				uint32 posCtx = (i * 8) / rec_.seqLen;
				q = dec_.illu8Coder->coder.DecodeSymbol(dec_.illu8Coder->rc, posCtx);
				ASSERT(q < 8);
			}

			rec_.qua[ii] = params_.archType.qualityOffset + idxToQua_8bin[q];
		}
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		struct cond_quantizer_list_t *qlist = globalQuaData.codebook.qlist;
		struct quantizer_t *q;
		uint32_t idx = 0;
		uint32_t qv_hat_prev = 0;
		uint32_t qv_hat = 0, q_state = 0;

		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 ii = rec_.IsReadReverse() ? rec_.seqLen - 1 - i : i;

			// decode the state [we prefer to compress the state of the quantizer than the actual value of qv_hat]
			// choose the quantizer given the current pos and the previous quantized value
			q = choose_quantizer(qlist, &local_well, i, qv_hat_prev, &idx);

			//qv_hat = dec_.qvzCoder->coder.DecodeSymbol(dec_.qvzCoder->rc);
			q_state = dec_.qvzCoder->DecodeNext(i, idx);
			qv_hat = q->output_alphabet->symbols[q_state];

			// get the value from state
			ASSERT(qv_hat < 64);

			// store the current qv_hat for next iteration
			qv_hat_prev = qv_hat;

			// reconstruct the original scale
			rec_.qua[ii] = qv_hat + params_.archType.qualityOffset;
		}
		break;
	}

	}
}


void IHeaderStoreBase::SetupFieldCompressionSpec()
{
	// select compression method(s) -- HINT: can be done globally per archive
	// or locally per block to further improve the ratio
	if (!fieldsSpec.empty())
		fieldsSpec.clear();
	fieldsSpec.resize(headData.fields.size());

	for (uint32 i = 0; i < headData.fields.size(); ++i)
	{
		const auto& f = headData.fields[i];
		auto& fcs = fieldsSpec[i];

		if (f.isConst)
		{
			fcs.method = FieldCompressionSpec::COMP_CONST;
		}
		else if (f.isNumeric)
		{
			fcs.method = FieldCompressionSpec::COMP_RAW;
			uint64 diff = f.maxValue - f.minValue;
			fcs.bitsPerValue = int_log(diff, 2) + 1;
		}
		else
		{
			fcs.method = FieldCompressionSpec::COMP_TOKEN;
			fcs.bitsPerValue = int_log(f.possibleValues.size(), 2) + 1;
			fcs.tokenValues.insert(fcs.tokenValues.end(), f.possibleValues.begin(), f.possibleValues.end());
		}
	}
}

void IHeaderStoreBase::CompressReadId(const FastqRecord& rec_, FieldEncoders& enc_) const
{
	uint32 fieldStartPos = 0;

	auto curFieldInfo = headData.fields.begin();
	auto curFieldComp = fieldsSpec.begin();

	for (uint32 i = 0; i <= rec_.headLen; ++i)
	{
		// TODO: go field by field, separator by separator
		if (rec_.head[i] != curFieldInfo->separator && (i != rec_.headLen))
			continue;

		// TODO: skip the field compression in PE case
		//

		if (curFieldComp->method == FieldCompressionSpec::COMP_CONST)
		{
			fieldStartPos = i + 1;
			curFieldComp++;
			curFieldInfo++;
			continue;
		}

		const uint32 fieldId = curFieldInfo - headData.fields.begin();
		const char* field = rec_.head + fieldStartPos;
		const uint32 fieldLen = i - fieldStartPos;

		uint64 val = 0;
		bool isNum = is_num(field, fieldLen, val);

		ASSERT((isNum && curFieldComp->method == FieldCompressionSpec::COMP_RAW)
			   || (!isNum && curFieldComp->method == FieldCompressionSpec::COMP_TOKEN));

		if (curFieldComp->method == FieldCompressionSpec::COMP_TOKEN)
		{
			auto it = std::find(curFieldComp->tokenValues.begin(), curFieldComp->tokenValues.end(), std::string(field, fieldLen));
			ASSERT(it != curFieldComp->tokenValues.end());

			uint32 id = it - curFieldComp->tokenValues.begin();
			enc_.tokenCoder->coder.EncodeSymbol(enc_.tokenCoder->rc, id, fieldId);
		}
		else
		{
			ASSERT(val >= curFieldInfo->minValue);

			int64 diff = val - curFieldInfo->minValue;

			uint32 ctxBase = (fieldId << 2);
			int32 valueRange = curFieldInfo->maxValue - curFieldInfo->minValue;

			// store the numeric values directly in arithmetic stream
			int32 plog = int_log(valueRange, 256);	// can be calculated once per field
			while (plog >= 0)
			{
				enc_.valueCoder->coder.EncodeSymbol(enc_.valueCoder->rc,
													(diff >> (8 * plog--)) & 0xFF,
													ctxBase++);
			}

			/*
			if (valueRange >= (1 << 24))
				enc_.valueCoder->coder.EncodeSymbol(enc_.valueCoder->rc, (diff >> 24) & 0xFF, ctxBase++);
			if (valueRange >= (1 << 16))
				enc_.valueCoder->coder.EncodeSymbol(enc_.valueCoder->rc, (diff >> 16) & 0xFF, ctxBase++);
			if (valueRange >= (1 << 8))
				enc_.valueCoder->coder.EncodeSymbol(enc_.valueCoder->rc, (diff >> 8) & 0xFF, ctxBase++);

			enc_.valueCoder->coder.EncodeSymbol(enc_.valueCoder->rc, diff & 0xFF, ctxBase);
			*/
		}


		fieldStartPos = i + 1;
		curFieldComp++;
		curFieldInfo++;
	}

	ASSERT(curFieldInfo == headData.fields.end());
}

void IHeaderStoreBase::DecompressReadId(FastqRecord& rec_, FieldDecoders& dec_) const
{
	// WARN: we won't know the header length in advance!

	auto curFieldInfo = headData.fields.begin();
	auto curFieldComp = fieldsSpec.begin();

	rec_.headLen = 0;

	for (uint32 i = 0; i < headData.fields.size(); ++i)
	{
		const uint32 fieldId = curFieldInfo - headData.fields.begin();

		switch (curFieldComp->method)
		{
		case FieldCompressionSpec::COMP_CONST:
		{
			if (curFieldInfo->isNumeric)
			{
				// copy the constant numeric value (min/max)
				uint32 len = to_string(rec_.head + rec_.headLen, curFieldInfo->minValue);
				ASSERT(len > 0);
				rec_.headLen += len;
			}
			else
			{
				// copy the only one unique value
				const auto& val = *curFieldInfo->possibleValues.begin();
				std::copy(val.begin(), val.end(), rec_.head + rec_.headLen);
				rec_.headLen += val.length();
			}

			break;
		}

		case FieldCompressionSpec::COMP_TOKEN:
		{
			uint32 id = dec_.tokenCoder->coder.DecodeSymbol(dec_.tokenCoder->rc, fieldId);
			ASSERT(id < curFieldComp->tokenValues.size());

			const auto& tok = curFieldComp->tokenValues[id];
			std::copy(tok.begin(), tok.end(), rec_.head + rec_.headLen);
			rec_.headLen += tok.length();

			break;
		}

		case FieldCompressionSpec::COMP_RAW:
		{
			int64 val = 0;
			int64 valueRange = curFieldInfo->maxValue - curFieldInfo->minValue;
			uint32 ctxBase = (fieldId << 2);

			while (valueRange > 0)
			{
				uint32 v = dec_.valueCoder->coder.DecodeSymbol(dec_.valueCoder->rc, ctxBase++);
				ASSERT(v < 256);

				val = (val << 8) | v;
				valueRange >>= 8;
			}

			/*
			if (valueRange >= (1 << 24))
				val = dec_.valueCoder->coder.DecodeSymbol(dec_.valueCoder->rc, ctxBase++) << 24;
			if (valueRange >= (1 << 16))
				val = val | (dec_.valueCoder->coder.DecodeSymbol(dec_.valueCoder->rc, ctxBase++) << 16);
			if (valueRange >= (1 << 8))
				val = val | (dec_.valueCoder->coder.DecodeSymbol(dec_.valueCoder->rc, ctxBase++) << 8);

			val = val | dec_.valueCoder->coder.DecodeSymbol(dec_.valueCoder->rc, ctxBase);
			*/

			val += curFieldInfo->minValue;
			ASSERT(val <= (int64)curFieldInfo->maxValue);

			uint32 len = to_string(rec_.head + rec_.headLen, val);
			ASSERT(len > 0);
			rec_.headLen += len;

			break;
		}
		}

		// do not add the separator if the field is the last one
		if (i != headData.fields.size() - 1)
		{
			rec_.head[rec_.headLen] = curFieldInfo->separator;
			rec_.headLen += 1;
		}

		curFieldInfo++;
		curFieldComp++;
	}

	ASSERT(rec_.headLen > 0);
}


void ILzCompressorBase::StoreHeader(const LzBlockHeader& header_, DataChunk &buffer_)
{
	BitMemoryWriter blockWriter(buffer_.data);

	StoreRawHeader(header_, blockWriter);

	// save work buffer sizes
	//
	for (uint64 s : header_.workBufferSizes)
		blockWriter.Put8Bytes(s);		// can reduce to 4B

	// save PPMd buffer sizes
	//
	for (uint64 s : header_.compBufferSizes)
		blockWriter.Put8Bytes(s);		// can reduce to 4B
}


void ILzCompressorBase::ReadHeader(LzBlockHeader& header_, DataChunk& buffer_)
{
	BitMemoryReader blockReader(buffer_.data, buffer_.size);

	ReadRawHeader(header_, blockReader);

	// HINT: this can be read in one loop -- we need to change only store order
	//
	uint64 totalBlockSize = LzBlockHeader::Size(header_.compBufferSizes.size());
	for (uint64& s : header_.workBufferSizes)
		s = blockReader.Get8Bytes();

	for (uint64& s : header_.compBufferSizes)
	{
		s = blockReader.Get8Bytes();
		totalBlockSize += s;
	}
	ASSERT(totalBlockSize + header_.footerSize == buffer_.size);
	if (header_.rawIdStreamSize != 0)
	{
		ASSERT(blockReader.Position() == LzBlockHeader::Size(header_.compBufferSizes.size()));
	}
	else
	{
		ASSERT(blockReader.Position() == LzBlockHeader::Size(header_.compBufferSizes.size()) - sizeof(header_.rawIdStreamSize));
	}
}


void IDnaRawStoreBase::StoreHeader(const RawBlockHeader &header_, DataChunk &buffer_)
{
	BitMemoryWriter blockWriter(buffer_.data);

	StoreRawHeader(header_, blockWriter);

	blockWriter.Put8Bytes(header_.dnaCompSize);
	blockWriter.Put8Bytes(header_.quaCompSize);

	if (params.archType.readsHaveHeaders)
	{
		blockWriter.Put8Bytes(header_.idTokenCompSize);
		blockWriter.Put8Bytes(header_.idValueCompSize);
	}
}


void IDnaRawStoreBase::ReadHeader(RawBlockHeader &header_, DataChunk &buffer_)
{
	BitMemoryReader blockReader(buffer_.data, buffer_.size);

	ReadRawHeader(header_, blockReader);

	header_.dnaCompSize = blockReader.Get8Bytes();
	header_.quaCompSize = blockReader.Get8Bytes();

	if (params.archType.readsHaveHeaders)
	{
		header_.idTokenCompSize = blockReader.Get8Bytes();
		header_.idValueCompSize = blockReader.Get8Bytes();
	}
	else
	{
		header_.idTokenCompSize = 0;
		header_.idValueCompSize = 0;
	}
}


// LzCompressorSE / LzDecompressorSE
//
//
LzCompressorSE::LzCompressorSE(const CompressorParams& params_,
							   const QualityCompressionData& globalQuaData_,
							   const FastqRawBlockStats::HeaderStats& headData_,
							   const CompressorAuxParams& auxParams_)
	:	ILzCompressorBase(params_, globalQuaData_, headData_, auxParams_)
	,	readsClassifier(params.minimizer, params_.classifier)
	,   blockStats(NULL)
	,	ppmdCoder(NULL)
	,	consBuilder(params_.consensus, params_.minimizer)
{
	ppmdCoder = new PpmdEncoder();
	ppmdCoder->StartCompress(DefaultPpmdOrder, DefaultPpmdMemorySizeMb);
}


LzCompressorSE::~LzCompressorSE()
{
	TFree(mainCtx.dna.lettersXRcCoder);
	TFree(mainCtx.dna.lettersCRcCoder);
	TFree(mainCtx.dna.revRcCoder);
	TFree(mainCtx.dna.matchRcCoder);
	TFree(mainCtx.dna.matchRleCoder);
	TFree(mainCtx.dna.lzRle0Coder);
#if (ENC_HR_AC)
	TFree(mainCtx.dna.hardCoder);
#endif

	ppmdCoder->FinishCompress();
	delete ppmdCoder;
}


void LzCompressorSE::SetupBuffers(Buffer& buffer_, std::vector<bool>& ppmdBufferCompMask_, uint32 bufferNum_)
{
	uint64 rawOutSize = std::accumulate(blockDesc.header.workBufferSizes.begin(),
										blockDesc.header.workBufferSizes.end(), 0);

	// usually the output data is compressed to 1/3-1/4 of the input data, but we want to keep the margin
	//
	uint64 preallocSize = rawOutSize;

	if (preallocSize > (1 << 25))
		preallocSize /= 2;
	else if (preallocSize <= (1 << 10))
		preallocSize *= 2;

	if (buffer_.Size() < preallocSize)
		buffer_.Extend(preallocSize);


	// setup buffer compression mask -- TODO: create only once and reuse
	//
	SetupBufferMask(ppmdBufferCompMask_, bufferNum_);
}


void LzCompressorSE::StartEncoding(std::vector<DataChunk*>& buffers_)
{
	ASSERT(mainCtx.writers.size() == 0);

	// create the main lz context
	//
	lzContextStack.push(LzContext());
	LzContext& lzCtx = lzContextStack.top();
	params.minimizer.GenerateMinimizer(blockDesc.header.minimizerId, lzCtx.currentMinimizerBuf.data());
	lzCtx.currentMinimizerBuf[params.minimizer.signatureLen] = 0;


	//  TODO: create only once in ctor() and use Reset()
	//

	// create mainCtx.writers
	//
	for (DataChunk* buf : buffers_)
	{
		mainCtx.writers.push_back(new BitMemoryWriter(buf->data));
		buf->size = 0;
	}

	// create dna encoders
	//
//	flagCoder = new FlagEncoder(*mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]);
	mainCtx.dna.revRcCoder = new DnaEncoders::RevEncoder(*mainCtx.writers[FastqWorkBuffersSE::RevBuffer]);
	mainCtx.dna.matchRcCoder = new DnaEncoders::RevEncoder(*mainCtx.writers[FastqWorkBuffersSE::MatchBinaryBuffer]);
	mainCtx.dna.lettersXRcCoder = new DnaEncoders::LettersEncoder(*mainCtx.writers[FastqWorkBuffersSE::LetterXBuffer]);
	mainCtx.dna.lettersCRcCoder = new DnaEncoders::LettersEncoder(*mainCtx.writers[FastqWorkBuffersSE::ConsensusLetterBuffer]);
	mainCtx.dna.matchRleCoder = new BinaryRleEncoder(*mainCtx.writers[FastqWorkBuffersSE::MatchBuffer]);
	mainCtx.dna.consMatchCoder = new BinaryRleEncoder(*mainCtx.writers[FastqWorkBuffersSE::ConsensusMatchBuffer]);
	mainCtx.dna.lzRle0Coder = new Rle0Encoder(*mainCtx.writers[FastqWorkBuffersSE::LzIdBuffer]);
#if (ENC_HR_AC)
	mainCtx.dna.hardCoder = new DnaEncoders::HardEncoder(*mainCtx.writers[FastqWorkBuffersSE::HardReadsBuffer]);
#endif

	// start encoders
	//
	mainCtx.dna.matchRleCoder->Start();
	mainCtx.dna.consMatchCoder->Start();
	mainCtx.dna.lzRle0Coder->Start();
//	flagCoder->Start();
	mainCtx.dna.revRcCoder->Start();
	mainCtx.dna.matchRcCoder->Start();
	mainCtx.dna.lettersXRcCoder->Start();
	mainCtx.dna.lettersCRcCoder->Start();
#if (ENC_HR_AC)
	mainCtx.dna.hardCoder->Start();
#endif

	// create/set quality encoders
	//
	{
		BitMemoryWriter* writer = mainCtx.writers[FastqWorkBuffersSE::QualityBuffer];
		switch (params.quality.method)
		{
		case QualityCompressionParams::MET_NONE:
		{
			mainCtx.qua.rawCoder = writer;	// just link the coder
			break;
		}

		case QualityCompressionParams::MET_BINARY:
		{
			mainCtx.qua.binaryCoder = new QualityEncoders::BinaryEncoder(*writer);
			mainCtx.qua.binaryCoder->Start();
			break;
		}

		case QualityCompressionParams::MET_8BIN:
		{
			mainCtx.qua.illu8Coder = new QualityEncoders::Illu8Encoder(*writer);
			mainCtx.qua.illu8Coder->Start();
			break;
		}

		case QualityCompressionParams::MET_QVZ:
		{
			mainCtx.qua.qvzCoder = new /*QualityEncoders::*/QVZEncoder(writer, globalQuaData.codebook.qlist);
			mainCtx.qua.qvzCoder->Start();

			// also reset rng when starting encoding so that the random number
			// generator will have the same initial seed also when
			// compressing in multithreaded mode
			ResetWellRng();
			break;
		}
		}
	}


	// create headers encoders
	//
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.id.tokenCoder = new FieldEncoders::TokenEncoder(*mainCtx.writers[FastqWorkBuffersSE::ReadIdTokenBuffer]);
		mainCtx.id.valueCoder = new FieldEncoders::ValueEncoder(*mainCtx.writers[FastqWorkBuffersSE::ReadIdValueBuffer]);

		mainCtx.id.tokenCoder->Start();
		mainCtx.id.valueCoder->Start();
	}


	// prepare for dry run
	//
	if (auxParams.dry_run)
	{
		dryFastqWriter->SetPosition(0);

		currentRecordIdx = 1;
	}
}


void LzCompressorSE::EndEncoding(std::vector<DataChunk*>& buffers_)
{
	// end dna encoders -- TODO: move this resposibility to DnaEncoders
	//
	mainCtx.dna.matchRleCoder->End();
	mainCtx.dna.consMatchCoder->End();
	mainCtx.dna.lzRle0Coder->End();
//	flagCoder->End();
	mainCtx.dna.revRcCoder->End();
	mainCtx.dna.matchRcCoder->End();
	mainCtx.dna.lettersXRcCoder->End();
	mainCtx.dna.lettersCRcCoder->End();
#if (ENC_HR_AC)
	mainCtx.dna.hardCoder->End();
#endif


	// end and free quality encoders -- TODO: move this responsibility to QualityEncoders
	//
	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		mainCtx.qua.rawCoder = NULL;		// this was just a link
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		mainCtx.qua.binaryCoder->End();
		TFree(mainCtx.qua.binaryCoder);
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		mainCtx.qua.illu8Coder->End();
		TFree(mainCtx.qua.illu8Coder);
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		mainCtx.qua.qvzCoder->End();
		TFree(mainCtx.qua.qvzCoder);

		break;
	}
	}


	// end read id encoders
	//
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.id.tokenCoder->End();
		mainCtx.id.valueCoder->End();

		TFree(mainCtx.id.tokenCoder);
		TFree(mainCtx.id.valueCoder);
	}


	// flush remaining data
	//
	for (uint32 i = 0; i < mainCtx.writers.size(); ++i)
	{
		auto writer = mainCtx.writers[i];
		writer->FlushPartialWordBuffer();
		buffers_[i]->size = writer->Position();
		blockDesc.header.workBufferSizes[i] = writer->Position();
	}

	// delete dna encoders -- TODO: create only once and reuse, move this resposnbility to DnaEncoders
	//
	TFree(mainCtx.dna.lettersXRcCoder);
	TFree(mainCtx.dna.lettersCRcCoder);
	TFree(mainCtx.dna.revRcCoder);
	TFree(mainCtx.dna.matchRcCoder);
//	TFree(flagCoder);
	TFree(mainCtx.dna.matchRleCoder);
	TFree(mainCtx.dna.consMatchCoder);
	TFree(mainCtx.dna.lzRle0Coder);
#if (ENC_HR_AC)
	TFree(mainCtx.dna.hardCoder);
#endif

	// delete writers
	//
	for (BitMemoryWriter* writer : mainCtx.writers)
		delete writer;

	mainCtx.writers.clear();

	lzContextStack.pop();
	ASSERT(lzContextStack.empty());


	// finalize for dry run
	//
	if (auxParams.dry_run)
	{
		dryFastqWriter->Flush();

		// POSIX requires C stdio FILE* ops to be atomic -- no mutex here required
		//
		fwrite(dryFastqWriter->Pointer(), 1, dryFastqWriter->Position(), auxParams.f_uncompressed);
	}
}


void LzCompressorSE::CompressBuffers(std::vector<DataChunk*>& buffers_,
										const std::vector<bool>& ppmdBufferCompMask_,
										DataChunk& compChunk_, uint64 chunkOffset_)
{
	byte* outMemBegin = compChunk_.data.Pointer();
	uint64 outMemPos = chunkOffset_;

	// copy range-coded buffers
	//
	for (uint32 i = 0; i < ppmdBufferCompMask_.size(); ++i)
	{
		if (ppmdBufferCompMask_[i])
			continue;

		DataChunk* buf = buffers_[i];
		if (buf->size > 0)
		{

			if (outMemPos + buf->size > compChunk_.data.Size())
			{
				compChunk_.data.Extend(outMemPos + buf->size, true);
				outMemBegin = compChunk_.data.Pointer();
			}

			std::copy(buf->data.Pointer(), buf->data.Pointer() + buf->size, outMemBegin + outMemPos);
			outMemPos += buf->size;
		}
		blockDesc.header.compBufferSizes[i] = buf->size;
	}


	// PPMd compress
	//
	for (uint32 i = 0; i < ppmdBufferCompMask_.size(); ++i)
	{
		if (!ppmdBufferCompMask_[i])
			continue;

		uint64 inSize = buffers_[i]->size;

		// take care of the case when we preallocated too small output memory
		if (inSize > 0)
		{
			// TODO: move to separate fuction
			//
			byte* inMem = buffers_[i]->data.Pointer();

			if (outMemPos + inSize > compChunk_.data.Size())
			{
				compChunk_.data.Extend(outMemPos + inSize, true);
				outMemBegin = compChunk_.data.Pointer();
			}

			byte* outMem = outMemBegin + outMemPos;
			uint64_t outSize = compChunk_.data.Size() - outMemPos;

			bool r = ppmdCoder->EncodeNextMember(inMem, inSize, outMem, outSize);
			ASSERT(r);
			ASSERT(outSize > 0);

			outMemPos += outSize;

			blockDesc.header.compBufferSizes[i] = outSize;
		}
		else
		{
			blockDesc.header.compBufferSizes[i] = 0;
		}
	}

	compChunk_.size = outMemPos;
}


void LzCompressorSE::Compress(const std::vector<FastqRecord>& reads_,
									PackContext& packCtx_,
									uint32 minimizerId_,
									uint64 rawDnaStreamSize_,
									FastqCompressedBin& fastqWorkBin_,
									CompressedFastqBlock &compBin_)
{
	ASSERT(reads_.size() > 0);
	ASSERT(packCtx_.graph->nodes.size() > 0);

	ASSERT(minimizerId_ != params.minimizer.SignatureN());
	ASSERT(packCtx_.stats.maxSeqLen > 0);
	ASSERT(packCtx_.stats.maxSeqLen >= packCtx_.stats.minSeqLen);


	// TODO: better handling of stats and log
	//
	blockStats = &compBin_.stats;
	blockStats->currentSignature = minimizerId_;
	//
	// //


	// initialize block description header
	//
	const uint32 buffersNum = fastqWorkBin_.buffers.size();
	blockDesc.Reset(buffersNum);
	blockDesc.header.minimizerId = minimizerId_;
	blockDesc.header.recordsCount = reads_.size();

	blockDesc.header.recMinLen = packCtx_.stats.minSeqLen;
	blockDesc.header.recMaxLen = packCtx_.stats.maxSeqLen;
	blockDesc.header.rawDnaStreamSize = rawDnaStreamSize_;
	blockDesc.header.rawIdStreamSize = 0;					// calculated no-fly; TODO: raw dna should be too



	// prepare mainCtx.writers and start encoding
	//
	StartEncoding(fastqWorkBin_.buffers);



	// TODO: shall we compress small and large bins in a different way???
	CompressRecords(packCtx_);
	ASSERT(blockStats->recordsCount == reads_.size() ||
		   printf("WARN: (%d) rc: %ld ; raw: %ld \n", minimizerId_, blockStats->recordsCount, reads_.size()) == 0);


	EndEncoding(fastqWorkBin_.buffers);


	// compress quality -- if additional processing is needed after encoding
	//
	CompressQuality();



	// setup and compress buffers
	//
	SetupBuffers(compBin_.dataBuffer.data, mainCtx.bufferCompMask, buffersNum);

	uint64 offset = LzBlockHeader::Size(buffersNum);
	CompressBuffers(fastqWorkBin_.buffers, mainCtx.bufferCompMask, compBin_.dataBuffer, offset);




	// store footer and header
	//
	compBin_.signatureId = minimizerId_;
	{
		blockDesc.header.footerOffset = compBin_.dataBuffer.size;

		BitMemoryWriter writer(compBin_.dataBuffer.data, compBin_.dataBuffer.size);
		StoreRawFooter(blockDesc.footer, writer);

		blockDesc.header.footerSize = writer.Position() - blockDesc.header.footerOffset;
		compBin_.dataBuffer.size = writer.Position();
	}

	StoreHeader(blockDesc.header, compBin_.dataBuffer);



	// save stats -- only for debug purposes
	//
	auto& cBuffer = compBin_.stats.bufferSizes["CompSize"];
	cBuffer.resize(blockDesc.header.compBufferSizes.size());
	for (uint32 i = 0; i < blockDesc.header.compBufferSizes.size(); ++i)
		cBuffer[i] = blockDesc.header.compBufferSizes[i];

	auto& rBuffer = compBin_.stats.bufferSizes["RawSize"];
	rBuffer.resize(blockDesc.header.compBufferSizes.size());
	for (uint32 i = 0; i < blockDesc.header.workBufferSizes.size(); ++i)
		rBuffer[i] = blockDesc.header.workBufferSizes[i];
}


void LzCompressorSE::CompressRecords(PackContext& packCtx_)
{
	// sort records
	//
	TFastqComparator<const MatchNode&> comparator;
	std::sort(packCtx_.graph->nodes.begin(), packCtx_.graph->nodes.end(), comparator);


	// match records
	//
	readsClassifier.ConstructMatchTree(*packCtx_.graph, packCtx_.rootNodes);


	// encode matching trees
	//
	for (MatchNode* root : packCtx_.rootNodes)
	{
		// try to build contigs traversing the tree
		//
		if (root->HasChildren())
			consBuilder.Build(root, packCtx_);


		// encode the tree
		//
		std::deque<MatchNode*> nq(1, root);
		while (!nq.empty())
		{
			MatchNode* node = nq.front();
			nq.pop_front();

			CompressNode(node);

			if (node->HasChildren())
				nq.insert(nq.end(), node->children->begin(), node->children->end());

			// we can clear the node after compression
			//
			node->Clear();
		}

		//uint32 ts = Consensus::TreeSize(newRoot, false);
		//blockStats->treeSizeFreqs[ts]++;

		// we can now clear the contigs
		//
		packCtx_.Clear(true);
	}
}


void LzCompressorSE::CompressNode(MatchNode* node_)
{
	LzContext& lzEncoder = lzContextStack.top();

	if (node_->type == MatchNode::TYPE_HARD)		// hard read - the tree root node
	{
		lzEncoder.history.clear();

#if EXTRA_MEM_OPT
		lzEncoder.history.shrink_to_fit();
#endif

		FastqRecord* rec = node_->record;

		CompressRead(*rec, HardRead);

		lzEncoder.history.push_back(rec);

		// compress exact matches and sub-trees (if present)
		if (node_->HasExactMatches())
			CompressExactChildren(node_);

		if (node_->HasSubTreeGroup())
			CompressSubTree(node_);


		// TODO: here we can already clear the groups to free some space
		//


		/*
		// update stats
		//
		const bool isNonEncoding = !node_->HasChildren();

		if (isNonEncoding)
			blockStats->orphanHardReads++;

		const bool encodesTree = node_->HasSubTreeGroup();
		if (encodesTree)
			blockStats->codingSubHardReads++;

		uint32 childCount = 0;
		if (node_->HasChildren())
			childCount = node_->children->size();

		blockStats->nodeStdChildrenFreqs[childCount]++;
		*/
	}
	else
	{
		// standard lz-match
		CompressMatch(node_, lzEncoder);


		// encode exact matches and sub trees (if present)
		//
		if (node_->HasExactMatches())
			CompressExactChildren(node_);

		if (node_->HasSubTreeGroup())
			CompressSubTree(node_);

		// here we can also already delete the groups
		//

		if (node_->HasContigGroup())	// consensus node
		{
			CompressContig(node_);

			// TODO: here we can already clear the groups to free some space
			//

		}


		/*
		// update stats
		//
		if (!node_->HasChildren())
			blockStats->orphanLzMatches++;

		if (node_->HasSubTreeGroup())
			blockStats->codingSubLzMatches++;

		uint32 childCount = 0;
		if (node_->HasChildren())
			childCount = node_->children->size();

		blockStats->nodeStdChildrenFreqs[childCount]++;
		//blockStats->mismFreqs[mismCount]++;
		blockStats->lzMatches++;
		*/
	}
}


void LzCompressorSE::CompressExactChildren(const MatchNode *node_)
{
	ASSERT(node_->HasExactMatches());

	const auto& emGroup = node_->GetExactMatches()->records;
	for (auto rec : emGroup)
	{
		CompressRead(*rec, ExactMatch);
	}
}


void LzCompressorSE::CompressExactRead(const FastqRecord &record_)
{
	mainCtx.dna.revRcCoder->coder.EncodeSymbol(mainCtx.dna.revRcCoder->rc, (uint32)record_.IsReadReverse());		// we lose this information when encoding CONS

	mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(ReadIdentical);


	// update stats
	//
	blockStats->counts["ExactMatches"]++;
}


void LzCompressorSE::CompressHardRead(const FastqRecord& record_)
{
	mainCtx.dna.revRcCoder->coder.EncodeSymbol(mainCtx.dna.revRcCoder->rc, (uint32)record_.IsReadReverse());

	mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(ReadDifficult);

	// CHECK ME
	//mainCtx.writers[FastqWorkBuffersSE::HardReadsBuffer]->PutBytes((byte*)rec.seq, rec.minimPos);
	//mainCtx.writers[FastqWorkBuffersSE::HardReadsBuffer]->PutByte(MinimizerPositionSymbol);
	//mainCtx.writers[FastqWorkBuffersSE::HardReadsBuffer]->PutBytes((byte*)rec.seq + rec.minimPos + params.minimizer.signatureLen, rec.seqLen);

	for (int32 i = 0; i < record_.seqLen; ++i)
	{
#if (ENC_HR_AC)
		if (i < record_.minimPos || i >= record_.minimPos + params.minimizer.signatureLen)
		{
			uint32 x = dnaToIdx[record_.seq[i]];
			ASSERT(x <= 4);
			mainCtx.dna.hardCoder->coder.EncodeSymbol(mainCtx.dna.hardCoder->rc, x);
		}
		else if (i == record_.minimPos)
			mainCtx.dna.hardCoder->coder.EncodeSymbol(mainCtx.dna.hardCoder->rc, 5);		// special symbol for m-mer
#else
		if (i < record_.minimPos || i >= record_.minimPos + params.minimizer.signatureLen)
		{
			mainCtx.writers[FastqWorkBuffersSE::HardReadsBuffer]->PutByte(record_.seq[i]);
		}
		else if (i == record_.minimPos)
		{
			mainCtx.writers[FastqWorkBuffersSE::HardReadsBuffer]->PutByte(MinimizerPositionSymbol);
		}
#endif
	}


	// update stats
	//
	blockStats->counts["HardReads"]++;
}

void LzCompressorSE::CompressNormalMatch(const MatchNode* node_,
											uint32 lzId_,
											bool expensiveEncoding_)
{
	const FastqRecord* rec = node_->record;
	const FastqRecord* lzRec = node_->lzRecord;

	// store rev-compl
	//
	mainCtx.dna.revRcCoder->coder.EncodeSymbol(mainCtx.dna.revRcCoder->rc, (uint32)rec->IsReadReverse());


	// store shift
	//
	ASSERT(ABS(node_->shiftValue) < 127);
	mainCtx.writers[FastqWorkBuffersSE::ShiftBuffer]->PutByte((int32)(ShiftOffset + node_->shiftValue));


	// store lz id
	//
	mainCtx.dna.lzRle0Coder->PutSymbol(lzId_);


	// store flag
	//
	int32 flag;
	if (node_->HasNoMismatches())
	{
		flag = ReadShiftOnly;
	}
	else
	{
		if (expensiveEncoding_)
			flag = ReadFullExpensive;
		else
			flag = ReadFullEncode;
	}
	mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(flag);



	// encode shift differences in respect to lz match
	//
	const char* bestSeq = lzRec->seq;
	uint32 bestLen = lzRec->seqLen;
	uint32 bestPos = lzRec->minimPos;

	const char* newSeq = rec->seq;
	uint32 newLen = rec->seqLen;

	if (node_->shiftValue >= 0)
	{
		ASSERT((int32)bestLen >= node_->shiftValue);
		ASSERT((int32)bestPos >= node_->shiftValue);

		bestSeq += node_->shiftValue;
		bestLen -= node_->shiftValue;
		bestPos -= node_->shiftValue;
	}
	else
	{
		ASSERT(newLen >= newLen);

		// encode insertion
		for (int32 i = 0; i < -node_->shiftValue; ++i)
		{
			mainCtx.dna.lettersXRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int32)rec->seq[i]], dnaToIdx[(int32)'N']);
		}

		newSeq += -node_->shiftValue;
		newLen -= -node_->shiftValue;

		// the '<=' was changed from '<' due to compression of rev-compl whole trees
		ASSERT((int32)bestPos - node_->shiftValue < (int32)(rec->seqLen));
	}
	uint32 minLen = MIN(bestLen, newLen);


	// encode differences in respect to lz match
	//
	if (flag == ReadFullEncode)
	{
		// TODO: optimize this routine
		for (uint32 i = 0; i < minLen; ++i)
		{
			if (i == bestPos)
			{
				i += params.minimizer.signatureLen - 1;
				continue;
			}

			if (bestSeq[i] == newSeq[i])
				mainCtx.dna.matchRleCoder->PutSymbol(true);
			else
			{
				mainCtx.dna.matchRleCoder->PutSymbol(false);
				mainCtx.dna.lettersXRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int32)newSeq[i]], dnaToIdx[(int32)bestSeq[i]]);
			}
		}
	}
	else if (flag == ReadFullExpensive)
	{
		for (uint32 i = 0; i < minLen; ++i)
		{
			if (i == bestPos)
			{
				i += params.minimizer.signatureLen - 1;
				continue;
			}

			mainCtx.dna.matchRcCoder->coder.EncodeSymbol(mainCtx.dna.matchRcCoder->rc, bestSeq[i] == newSeq[i]);

			if (bestSeq[i] != newSeq[i])
			{
				mainCtx.dna.lettersXRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int32)newSeq[i]], dnaToIdx[(int32)bestSeq[i]]);
			}
		}
	}


	// encode last shift insertion
	//
	for (uint32 i = minLen; i < newLen; ++i)
	{
		mainCtx.dna.lettersXRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int32)newSeq[i]], dnaToIdx[(int32)'N']);
	}


	// update stats
	//
	blockStats->counts["LzMatches"]++;
	blockStats->freqs["LzId"][lzId_]++;
	blockStats->freqs["Shift"][node_->shiftValue]++;
	blockStats->freqs["Cost"][node_->encodeCost]++;
}


void LzCompressorSE::CompressContig(MatchNode* node_)
{
	ASSERT(node_->HasContigGroup());


	// TODO: here we can already clear the groups to free some space
	const ContigDefinition* contig = node_->GetContigGroup();


	// prepare the consensus configuration
	//
	consEncoders.push(ConsensusEncoder());
	ConsensusEncoder& consEncoder = consEncoders.top();
	consEncoder.type = ConsensusEncoder::ConsensusGroup;
	consEncoder.definition = &contig->consensus;
	consEncoder.lastMinimPos = node_->record->minimPos;


	// start encoding
	//
	mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(ReadContigGroupStart);

	StoreContigDefinition(*consEncoder.definition,
						  node_->record->seq,
						  node_->record->minimPos);

	CompressContigRecords(node_, *contig);


	// end encoding
	//
	mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(ReadGroupEnd);


	// pop the consensus encoder
	//
	consEncoders.pop();


	// update stats
	//
	blockStats->counts["Contigs"]++;
	blockStats->freqs["ContigVariants"][contig->consensus.variantsCount]++;
	blockStats->freqs["ContigSizes"][contig->nodes.size()]++;

}


void LzCompressorSE::StoreContigDefinition(const ConsensusDefinition& definition_,
											  const char* mainSeq_,
											  uint32 mainSignaturePos_)
{
	// calculate the consensus lz range
	//
	std::pair<uint32, uint32> lzRefRange;						// we can calculate this in post-processing
	lzRefRange.first = definition_.readLen - mainSignaturePos_;
	lzRefRange.second = lzRefRange.first + definition_.readLen;

	if (definition_.readLen < 128)
	{
		const int32 r1Rescale = ShiftOffset - definition_.readLen / 2;
		const int32 r2Rescale = ShiftOffset - definition_.readLen * 3 / 2;

		mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(definition_.range.first + r1Rescale);
		mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(definition_.range.second + r2Rescale);
	}
	else
	{
		mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(definition_.range.first);

		// delta = range.second - range.first - seqLen + rescale
		int32 rangeDiff = (int32)definition_.range.second - (int32)definition_.range.first;
		int32 rescale = params.minimizer.signatureLen + ReadsContigBuilderParams::Default::BeginCut + ReadsContigBuilderParams::Default::EndCut;
		int32 delta = rangeDiff - (int32)definition_.readLen + rescale;

		ASSERT(delta > 0 && delta < 256);
		mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(delta);
	}


	// encode consensus according to last LZ match
	//
	for (uint32 i = definition_.range.first; i < definition_.range.second; ++i)
	{
		// skip the signature position range
		//
		if (i == definition_.readLen)
		{
			i += params.minimizer.signatureLen - 1;
			continue;
		}

		// encode the match mask
		//
		mainCtx.dna.consMatchCoder->PutSymbol(definition_.variantPositions[i] == 0);

		// store the missing consensus sequence
		//
		if (i < lzRefRange.first + ReadsContigBuilderParams::Default::BeginCut
				|| i >= lzRefRange.second - ReadsContigBuilderParams::Default::EndCut
				|| definition_.variantPositions[i] != 0)
		{
			mainCtx.dna.lettersCRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersCRcCoder->rc,
														dnaToIdx[(int)definition_.sequence[i]],
														dnaToIdx[(int)'N']);
		}
		else
		{
			ASSERT(definition_.sequence[i] == mainSeq_[i - lzRefRange.first]);
		}
	}
}


void LzCompressorSE::CompressContigRecords(MatchNode* /*mainNode_*/,
										   const ContigDefinition& contig_)
{
	// micro-optimization -- first record shift saved into other shift buffer
	bool firstRecord = true;


	LzContext& lzEncoder = lzContextStack.top();

	for (const MatchNode* match : contig_.nodes)
	{
		if (!firstRecord)
			mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(ReadContigGroupNext);

		CompressRead(*match->record, ContigRead, firstRecord);
		firstRecord = false;


		// update the LZ search window and stats
		//
		lzEncoder.history.push_back(match->record);


		// store the EMs or subtrees (if present)
		//
		if (match->HasExactMatches())
			CompressExactChildren(match);

		if (match->HasSubTreeGroup())
			CompressSubTree(match);
	}

	// here we could already free the consensus taking care of lz-matches
	// i.e. freeing the non-encoding ones
}


void LzCompressorSE::CompressSubTree(const MatchNode* node_)
{
	ASSERT(node_->HasSubTreeGroup());

	const FastqRecord& rootRec = *node_->record;
	const auto& treeList = node_->GetSubTrees();

	for (GraphEncodingContext* tree : treeList)
	{
		consEncoders.push(ConsensusEncoder());
		ConsensusEncoder& consEncoder = consEncoders.top();

		consEncoder.type = ConsensusEncoder::TreeGroup;
		consEncoder.lastMinimPos = tree->mainSignaturePos;


#if DEV_DEBUG_MODE
		// DEBUG: just check the signature correctness
		//
		std::array<char, MAX_SIGNATURE_LEN> sigSeq;
		params.minimizer.GenerateMinimizer(tree->signatureId, sigSeq.data());

		ASSERT(std::search(rootRec.seq + consEncoder.lastMinimPos,
						   rootRec.seq + consEncoder.lastMinimPos + params.minimizer.signatureLen,
						   sigSeq.data(),
						   sigSeq.data() + params.minimizer.signatureLen) == rootRec.seq + consEncoder.lastMinimPos
				|| printf("WARN: mainRecord: %.100s \nSigSeq: %.8s \nlastMPos: %d\n",
						   rootRec.seq,
						   sigSeq.data(),
						   consEncoder.lastMinimPos) == 0);
#endif


		// store tree meta-data to synchronise the signature position
		//
		if (rootRec.seqLen * 2 >= 256)
		{
			// TODO: verify for variable-length reads
			ASSERT(rootRec.seqLen < 256);
			mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(consEncoder.lastMinimPos);
		}
		else
		{
			int32 shift = consEncoder.lastMinimPos - (int32)rootRec.minimPos;
			mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(ShiftOffset + shift);
		}


		// start encoding
		//
		mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(ReadTreeGroupStart);


		// prepare lz context
		//
		lzContextStack.push(LzContext());
		LzContext& lzEncoder = lzContextStack.top();


		// prepare a local copy of the main read with altered minimizer position
		// which will serve as our new root
		//
		FastqRecord localMainRec = rootRec;
		localMainRec.minimPos = consEncoder.lastMinimPos;

		MatchNode localRootNode;
		localRootNode.record = &localMainRec;
		localRootNode.type = MatchNode::TYPE_NONE;		// important here for LZ-matching

		lzEncoder.history.push_back(&localMainRec);



		// create matching graph
		//
		/*
		if (params.useStoredTopology)
		{
			FastqReadsGraph::CreateGraphFromTopologyVector(graph.nodes,
														   tree.topology,
														   &localMainRec,
														   tree.records,
														   readsClassifier);
			graph.rootNodes.push_back(graph.nodes.front());
		}
		else
		*/
		PackContext localPackCtx(tree);

		readsClassifier.ConstructMatchTree(*tree, localPackCtx.rootNodes, &localRootNode);
		//ASSERT(localPackCtx.rootNodes.size() == 1);



		// iterate over root nodes... however here we should only have one root node, but
		// in some cases (eg re-bin with lower encode threshold) as a result of matching
		// one can obtain more
		//
		bool isFirstNode = true;

		for (MatchNode* root : localPackCtx.rootNodes)
		{
			//

			// try to build contigs
			//
			if (root->HasChildren())		// && tree size > min contig size !!!
				consBuilder.Build(root, localPackCtx);


			// in case of compressing the main node (first node), we only need to add its children
			// to the processing queue as the node itself has already been compressed
			MatchNodeIterator nodeIterator(root);
			std::deque<MatchNode*> nq;
			if (isFirstNode)
			{
				ASSERT(root->HasChildren());
				nq.insert(nq.end(), root->children->begin(), root->children->end());

				isFirstNode = false;
			}
			else
			{
				nq.push_back(root);
			}


			while (!nq.empty())
			{
				MatchNode* node = nq.front();
				nq.pop_front();

				CompressNode(node);

				if (node->HasChildren())
					nq.insert(nq.end(), node->children->begin(), node->children->end());


				// we can clear the node after compression
				//
				node->Clear();
			}


			// update sub-tree stats
			//uint32 ts = Consensus::TreeSize(newRoot, false);
			//blockStats->subTreeSizeFreqs[ts]++;
		}


		// end encoding
		//
		mainCtx.writers[FastqWorkBuffersSE::FlagBuffer]->PutByte(ReadGroupEnd);


		// clear the encoders
		//
		consEncoders.pop();
		lzContextStack.pop();


		// we can also clean the tree, as we won't need this data anymore
		//
		tree->Clear();
	}
}


void LzCompressorSE::CompressContigRead(const FastqRecord& rec_, bool useTreeShiftBuffer_)
{
	// delta encode shifts according to previous read
	//
	ConsensusEncoder& consEncoder = consEncoders.top();		// TODO: pass by arguments

	int32 dpos = (int32)rec_.minimPos - (int32)consEncoder.lastMinimPos;
	consEncoder.lastMinimPos = (int32)rec_.minimPos;

#if 0
	int32 delta = ShiftOffset + dpos;
	if (!useTreeShiftBuffer_)
	{
		ASSERT(delta < 256);
		mainCtx.writers[FastqWorkBuffersSE::ConsensusShiftBuffer]->PutByte(delta);
	}
	else
	{
		if (rec_.seqLen * 2 >= 256)
		{	
			ASSERT(rec_.seqLen < 256);
			mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(rec_.minimPos);

		}
		else
		{
			ASSERT(delta < 256);
			mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer]->PutByte(delta);
		}
	}
#else
	BitMemoryWriter* writer = NULL;

	if (!useTreeShiftBuffer_)
		writer = mainCtx.writers[FastqWorkBuffersSE::ConsensusShiftBuffer];
	else
		writer = mainCtx.writers[FastqWorkBuffersSE::TreeShiftBuffer];

	int32 delta = ShiftOffset + dpos;
	if (rec_.seqLen * 2 >= 256)
	{
		ASSERT(rec_.seqLen < 256);
		writer->PutByte(rec_.minimPos);
	}
	else
	{
		ASSERT(delta < 256);
		writer->PutByte(delta);
	}

#endif

	// store rev-compl info
	//
	mainCtx.dna.revRcCoder->coder.EncodeSymbol(mainCtx.dna.revRcCoder->rc, (uint32)rec_.IsReadReverse());


	// calculate the read range inside the consensus
	//
	const uint32 readLen = rec_.seqLen;
	const uint32 consStart = readLen - rec_.minimPos;

	ASSERT(consStart + ReadsContigBuilderParams::Default::BeginCut >= consEncoder.definition->range.first
		   || (rec_.minimPos <= ReadsContigBuilderParams::Default::BeginCut
			   && consStart + rec_.minimPos + params.minimizer.signatureLen >= consEncoder.definition->range.first));
	ASSERT(consStart + readLen - ReadsContigBuilderParams::Default::EndCut <= consEncoder.definition->range.second
		   || ((int32)rec_.seqLen - rec_.minimPos - params.minimizer.signatureLen <= (int32)ReadsContigBuilderParams::Default::EndCut
			   && consStart + rec_.minimPos - ReadsContigBuilderParams::Default::EndCut <= consEncoder.definition->range.second)
		   || printf("WARN: min pos: %d ; cons start: %d ; cons end: %d \n",
									  rec_.minimPos,
									  consEncoder.definition->range.first,
									  consEncoder.definition->range.second) == 0);


	// store the begin cut of the read
	//
	uint32 encodeIterator = 0;
	while (encodeIterator < ReadsContigBuilderParams::Default::BeginCut)
	{
		if (encodeIterator == rec_.minimPos)
		{
			ASSERT(encodeIterator + params.minimizer.signatureLen > ReadsContigBuilderParams::Default::BeginCut);
			encodeIterator += params.minimizer.signatureLen;
			continue;
		}

		mainCtx.dna.lettersCRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersCRcCoder->rc,
										   dnaToIdx[(int)rec_.seq[encodeIterator]],
										   dnaToIdx[(int)consEncoder.definition->sequence[consStart + encodeIterator]]);
		encodeIterator++;
	}


	// store the main part of the read
	//
	while (encodeIterator < readLen - ReadsContigBuilderParams::Default::EndCut)
	{
		// skip encoding of the information inside the signature
		//
		if (encodeIterator == rec_.minimPos)
		{
			encodeIterator += params.minimizer.signatureLen;
			continue;
		}

		if (consEncoder.definition->variantPositions[consStart + encodeIterator])
		{
			mainCtx.dna.lettersCRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersCRcCoder->rc,
											   dnaToIdx[(int)rec_.seq[encodeIterator]],
											   dnaToIdx[(int)consEncoder.definition->sequence[consStart + encodeIterator]]);
		}
		else
		{
			ASSERT(rec_.seq[encodeIterator] == consEncoder.definition->sequence[consStart + encodeIterator]);
		}
		encodeIterator++;
	}


	// store the end cut
	//
	while (encodeIterator < readLen)
	{
		mainCtx.dna.lettersCRcCoder->coder.EncodeSymbol(mainCtx.dna.lettersCRcCoder->rc,
										   dnaToIdx[(int)rec_.seq[encodeIterator]],
										   dnaToIdx[(int)consEncoder.definition->sequence[consStart + encodeIterator]]);
		encodeIterator++;
	}


	// update stats
	//
	blockStats->counts["ContigReads"]++;
}



void LzCompressorSE::CompressMatch(const MatchNode* node_, LzContext& lzEncoder_)
{
	if (auxParams.dry_run)
	{
		DryRunCompress(*node_->record, dryFastqWriter);
	}
	else
	{
		// compress node
		//
		ASSERT(node_->type == MatchNode::TYPE_LZ);

		uint32 lzParentIdx = std::find(lzEncoder_.history.rbegin(),
									   lzEncoder_.history.rend(),
									   node_->lzRecord) - lzEncoder_.history.rbegin();
		ASSERT(lzParentIdx != lzEncoder_.history.size());

		const uint32 mismCount = (node_->encodeCost - ABS(node_->shiftValue * params.classifier.shiftCost)) / params.classifier.mismatchCost;
		const bool isExpensiveEncoding = mismCount > params.maxMismatchesLowCost;

		CompressNormalMatch(node_, lzParentIdx, isExpensiveEncoding);

		lzEncoder_.history.push_back(node_->record);


		// compress header
		//
		if (params.archType.readsHaveHeaders)
		{
			CompressReadId(*node_->record, mainCtx.id);
			blockDesc.header.rawIdStreamSize += node_->record->headLen;
		}


		// compress quality
		//
		CompressReadQuality(*node_->record);


		blockStats->freqs["LzMatches-mism"][mismCount]++;
	}

	// update stats
	//
	blockStats->recordsCount++;
}


void LzCompressorSE::CompressRead(const FastqRecord& record_, ReadMatchType type_, bool aux_)
{
	if (auxParams.dry_run)
	{
		DryRunCompress(record_, dryFastqWriter);
	}
	else
	{
		// compress head data
		//
		if (params.archType.readsHaveHeaders)
		{
			CompressReadId(record_, mainCtx.id);
			blockDesc.header.rawIdStreamSize += record_.headLen;
		}


		// compress sequence data
		//
		switch (type_)
		{
		case HardRead: CompressHardRead(record_); break;
		case ExactMatch: CompressExactRead(record_); break;
		case ContigRead: CompressContigRead(record_, aux_); break;
		default:	ASSERT(0);
		}


		// compress quality data
		//
		CompressReadQuality(record_);
	}


	// update stats
	//
	blockStats->recordsCount++;
}

void LzCompressorSE::DryRunCompress(const FastqRecord &record_, BitMemoryWriter* dryFastqWriter_)
{
	ASSERT(auxParams.f_uncompressed != NULL);

	if (auxParams.output_fastq)
	{
		// output reads header
		//
		if (params.archType.readsHaveHeaders)
		{
			dryFastqWriter_->PutBytes((byte*)record_.head, record_.headLen);
		}
		else // dummy header -- this should not be the case
		{
			//const std::string libName = "SRX000000.X";

			char buf[64];
			int len = 0;
			buf[len++] = '@';
			params.minimizer.GenerateMinimizer(blockStats->currentSignature, buf + len);
			len += params.minimizer.signatureLen;
			buf[len++] = '.';
			len += to_string(buf + len, currentRecordIdx);

			dryFastqWriter_->PutBytes((byte*)buf, len);
		}
		dryFastqWriter_->PutByte('\n');


		// output dna
		//
		if (record_.IsReadReverse())
		{
			const char* codes = FastqRecord::GetRCCodes();
			for (int32 i = record_.seqLen - 1; i >= 0; i--)
				dryFastqWriter_->PutByte(codes[record_.seq[i] - 64]);
		}
		else
		{
			dryFastqWriter_->PutBytes((byte*)record_.seq, record_.seqLen);
		}
		dryFastqWriter_->PutByte('\n');

		// otput control line
		dryFastqWriter_->PutByte('+');
		dryFastqWriter_->PutByte('\n');

		currentRecordIdx++;
	}


	// output quality
	//
	{
		std::vector<byte> qscores;
		IQualityStoreBase::CompressReadQuality(record_, mainCtx.qua, params, true, &qscores);

		dryFastqWriter_->PutBytes(qscores.data(), qscores.size());
		dryFastqWriter_->PutByte('\n');
	}
}


void LzCompressorSE::CompressReadQuality(const FastqRecord& rec_)
{
	IQualityStoreBase::CompressReadQuality(rec_, mainCtx.qua, params);
}


void LzCompressorSE::CompressQuality()
{
	// HINT: here can be performed extra quality information compressed
	//
	// ...
	//
}



// decompressor
//
//
LzDecompressorSE::LzDecompressorSE(const CompressorParams &params_,
								   const QualityCompressionData& globalQuaData_,
								   const FastqRawBlockStats::HeaderStats& headData_)
	:	ILzCompressorBase(params_, globalQuaData_, headData_)
{
	// TODO: refactor -- we can initialize all mainCtx.readers and encoders in ctor
	//
	ppmdCoder = new PpmdDecoder();
	ppmdCoder->StartDecompress(DefaultPpmdMemorySizeMb);
}



LzDecompressorSE::~LzDecompressorSE()
{
	TFree(mainCtx.dna.matchRleCoder);
	TFree(mainCtx.dna.consMatchCoder);
	TFree(mainCtx.dna.lettersXRcCoder);
	TFree(mainCtx.dna.lettersCRcCoder);
	TFree(mainCtx.dna.revRcCoder);
	TFree(mainCtx.dna.matchRcCoder);
	//TFree(flagCoder);
	TFree(mainCtx.dna.lzRle0Coder);
#if (ENC_HR_AC)
	TFree(mainCtx.dna.hardCoder);
#endif

	ppmdCoder->FinishDecompress();
	delete ppmdCoder;
}


void LzDecompressorSE::StartDecoding(std::vector<DataChunk*>& buffers_)
{
	// prepare main lz context
	//
	lzContextStack.push(LzContext());
	LzContext& lzCtx = lzContextStack.top();
	params.minimizer.GenerateMinimizer(blockDesc.header.minimizerId, lzCtx.currentMinimizerBuf.data());
	lzCtx.currentMinimizerBuf[params.minimizer.signatureLen] = 0;


	// create mainCtx.readers -- TODO: create only once and reuse
	//
	ASSERT(mainCtx.readers.size() == 0);
	for (DataChunk* buf : buffers_)
	{
		mainCtx.readers.push_back(new BitMemoryReader(buf->data, buf->size));
	}


	// create decoders -- TODO: create only once and reuse
	//
	//flagCoder = new FlagDecoder(*mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]);
	mainCtx.dna.revRcCoder = new DnaDecoders::RevDecoder(*mainCtx.readers[FastqWorkBuffersSE::RevBuffer]);
	mainCtx.dna.matchRcCoder = new DnaDecoders::RevDecoder(*mainCtx.readers[FastqWorkBuffersSE::MatchBinaryBuffer]);
	mainCtx.dna.lettersXRcCoder = new DnaDecoders::LettersDecoder(*mainCtx.readers[FastqWorkBuffersSE::LetterXBuffer]);
	mainCtx.dna.lettersCRcCoder = new DnaDecoders::LettersDecoder(*mainCtx.readers[FastqWorkBuffersSE::ConsensusLetterBuffer]);
	mainCtx.dna.matchRleCoder = new BinaryRleDecoder(*mainCtx.readers[FastqWorkBuffersSE::MatchBuffer]);
	mainCtx.dna.consMatchCoder = new BinaryRleDecoder(*mainCtx.readers[FastqWorkBuffersSE::ConsensusMatchBuffer]);
	mainCtx.dna.lzRle0Coder = new Rle0Decoder(*mainCtx.readers[FastqWorkBuffersSE::LzIdBuffer]);
#if (ENC_HR_AC)
	mainCtx.dna.hardCoder = new DnaDecoders::HardDecoder(*mainCtx.readers[FastqWorkBuffersSE::HardReadsBuffer]);
#endif

	// start decoders
	//
	mainCtx.dna.matchRleCoder->Start();
	mainCtx.dna.consMatchCoder->Start();
	mainCtx.dna.lzRle0Coder->Start();
	//flagCoder->Start();
	mainCtx.dna.revRcCoder->Start();
	mainCtx.dna.matchRcCoder->Start();
	mainCtx.dna.lettersXRcCoder->Start();
	mainCtx.dna.lettersCRcCoder->Start();
#if (ENC_HR_AC)
	mainCtx.dna.hardCoder->Start();
#endif


	// create/set quality encoders
	//
	{
		BitMemoryReader* reader = mainCtx.readers[FastqWorkBuffersSE::QualityBuffer];
		switch (params.quality.method)
		{
		case QualityCompressionParams::MET_NONE:
		{
			mainCtx.qua.rawCoder = reader;	// just link the coder
			break;
		}

		case QualityCompressionParams::MET_BINARY:
		{
			mainCtx.qua.binaryCoder = new QualityDecoders::BinaryDecoder(*reader);
			mainCtx.qua.binaryCoder->Start();
			break;
		}

		case QualityCompressionParams::MET_8BIN:
		{
			mainCtx.qua.illu8Coder = new QualityDecoders::Illu8Decoder(*reader);
			mainCtx.qua.illu8Coder->Start();
			break;
		}

		case QualityCompressionParams::MET_QVZ:
		{	
			mainCtx.qua.qvzCoder = new /*QualityDecoders::*/QVZDecoder(reader, globalQuaData.codebook.qlist);
			mainCtx.qua.qvzCoder->Start();


			// also reset rng when starting decoding so that the random number
			// generator will have the same initial seed also when
			// compressing in multithreaded mode
			ResetWellRng();
			break;
		}
		}
	}

	// set fields encoder
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.id.tokenCoder = new FieldDecoders::TokenDecoder(*mainCtx.readers[FastqWorkBuffersSE::ReadIdTokenBuffer]);
		mainCtx.id.valueCoder = new FieldDecoders::ValueDecoder(*mainCtx.readers[FastqWorkBuffersSE::ReadIdValueBuffer]);

		mainCtx.id.tokenCoder->Start();
		mainCtx.id.valueCoder->Start();
	}
}


void LzDecompressorSE::EndDecoding()
{
	// end dna decoders
	//
	mainCtx.dna.matchRleCoder->End();
	mainCtx.dna.consMatchCoder->End();
	mainCtx.dna.lzRle0Coder->End();
	//flagCoder->End();
	mainCtx.dna.revRcCoder->End();
	mainCtx.dna.matchRcCoder->End();
	mainCtx.dna.lettersXRcCoder->End();
	mainCtx.dna.lettersCRcCoder->End();
#if (ENC_HR_AC)
	mainCtx.dna.hardCoder->End();
#endif


	// end quality decoders -- TODO: move this resposibility and reuse
	//
	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		mainCtx.qua.rawCoder = NULL;		// this was just a link
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		mainCtx.qua.binaryCoder->End();
		TFree(mainCtx.qua.binaryCoder);
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		mainCtx.qua.illu8Coder->End();
		TFree(mainCtx.qua.illu8Coder);
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		mainCtx.qua.qvzCoder->End();
		TFree(mainCtx.qua.qvzCoder);

		break;
	}
	}


	// finish fields encoder
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.id.tokenCoder->End();
		mainCtx.id.valueCoder->End();

		TFree(mainCtx.id.tokenCoder);
		TFree(mainCtx.id.valueCoder);
	}



	// delete decoders -- TODO: create only once and reuse
	//
	TFree(mainCtx.dna.matchRleCoder);
	TFree(mainCtx.dna.consMatchCoder);
	TFree(mainCtx.dna.lzRle0Coder);
	TFree(mainCtx.dna.lettersXRcCoder);
	TFree(mainCtx.dna.lettersCRcCoder);
	TFree(mainCtx.dna.revRcCoder);
	TFree(mainCtx.dna.matchRcCoder);
	//TFree(flagCoder);
#if (ENC_HR_AC)
	TFree(mainCtx.dna.hardCoder);
#endif


	for (uint32 i = 0; i < mainCtx.readers.size(); ++i)
		delete mainCtx.readers[i];
	mainCtx.readers.clear();


	// clear the main lz context
	//
	lzContextStack.pop();
	ASSERT(lzContextStack.empty());
}


void LzDecompressorSE::SetupBuffers(Buffer& outBuffer_, std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_)
{
	// setup output buffer
	//
	uint64 bufSize = blockDesc.header.rawDnaStreamSize * 2;		// * 2 -- for qualities
	if (params.archType.readsHaveHeaders)
		bufSize += blockDesc.header.rawIdStreamSize;

	if (outBuffer_.Size() < bufSize)
		outBuffer_.Extend(bufSize);

	SetupBufferMask(ppmdBufferCompMask_, buffersNum_);
}


void LzDecompressorSE::DecompressBuffers(std::vector<DataChunk*>& buffers_,
										   const std::vector<bool>& ppmdBufferCompMask_,
										   DataChunk &compChunk_,
										   uint64 chunkOffset_)
{
	// prepare buffers
	//
	for (uint32 i = 0; i < buffers_.size(); ++i)
	{
		if (buffers_[i]->data.Size() < blockDesc.header.workBufferSizes[i])
			buffers_[i]->data.Extend(blockDesc.header.workBufferSizes[i] + (blockDesc.header.workBufferSizes[i] / 8));
	}

	byte* inMemBegin = compChunk_.data.Pointer();
	uint64 inMemPos = chunkOffset_;


	// copy explicitly arithmetically encoded buffers
	//
	for (uint32 i = 0; i < ppmdBufferCompMask_.size(); ++i)
	{
		if (ppmdBufferCompMask_[i])
			continue;

		const uint64 size = blockDesc.header.workBufferSizes[i];
		DataChunk* buf = buffers_[i];

		if (size > 0)
		{
			if (buf->data.Size() < size)
				buf->data.Extend(size);

			std::copy(inMemBegin + inMemPos,
					  inMemBegin + inMemPos + size,
					  buf->data.Pointer());

			inMemPos += size;
		}

		buf->size = size;
	}


	// PPMd decompress
	//
	for (uint32 i = 0; i < ppmdBufferCompMask_.size(); ++i)
	{
		if (!ppmdBufferCompMask_[i])
			continue;

		uint64 inSize = blockDesc.header.compBufferSizes[i];

		if (inSize > 0)
		{
			// TODO: move to separate function
			//
			byte* inMem = inMemBegin + inMemPos;

			uint64_t outSize = buffers_[i]->data.Size();
			byte* outMem = buffers_[i]->data.Pointer();

			bool r = ppmdCoder->DecodeNextMember(inMem, inSize, outMem, outSize);
			ASSERT(r);
			ASSERT(outSize > 0);
			inMemPos += inSize;

			ASSERT(blockDesc.header.workBufferSizes[i] == outSize);
			buffers_[i]->size = outSize;
		}
		else
		{
			buffers_[i]->size = 0;
		}
	}
}


void LzDecompressorSE::Decompress(CompressedFastqBlock &compBin_,
										std::vector<FastqRecord>& reads_,
										FastqCompressedBin& fastqWorkBin_,
										FastqChunk& dnaBuffer_)
{
	// read header and footer
	//
	const uint32 buffersNum = fastqWorkBin_.buffers.size();
	blockDesc.Reset(buffersNum);
	ReadHeader(blockDesc.header, compBin_.dataBuffer);
	blockDesc.isLenConst = (blockDesc.header.recMinLen == blockDesc.header.recMaxLen);

	{
		BitMemoryReader reader(compBin_.dataBuffer.data,
							   blockDesc.header.footerOffset + blockDesc.header.footerSize,
							   blockDesc.header.footerOffset);
		ReadRawFooter(blockDesc.footer, reader);
	}



	// decompress quality if needed
	//
	DecompressQuality();


	// setup and decompress buffers
	//
	SetupBuffers(dnaBuffer_.data, mainCtx.bufferCompMask, buffersNum);


	DecompressBuffers(fastqWorkBin_.buffers, mainCtx.bufferCompMask,
					  compBin_.dataBuffer, LzBlockHeader::Size(buffersNum));

	// start decoding
	//
	StartDecoding(fastqWorkBin_.buffers);


	// prepare the output array and bind the buffers
	//
	reads_.resize(blockDesc.header.recordsCount);

	// -- decompress here the headers -- like in DSRC ;---)
	DecompressMetaAndLinkRecords(reads_, dnaBuffer_);

	//dnaBuffer_.size = blockDesc.header.rawDnaStreamSize * 2;


	// decode records
	//
	DecompressRecords(reads_);


	// finish decoding
	//
	EndDecoding();

	ASSERT(dnaBuffer_.size == blockDesc.header.rawDnaStreamSize * 2 + blockDesc.header.rawIdStreamSize);
}


void LzDecompressorSE::DecompressRecords(std::vector<FastqRecord>& reads_)
{
	uint64 recIdx = 0;
	while (recIdx < reads_.size())
	{
		ASSERT(consDecoders.size() == 0);

		uint32 flag = mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->GetByte();
		ASSERT(flag <= ReadContigGroupStart);

		// INFO: we are incrementing the read idx only when decompressing normal reads
		// the groups: sub-tree and contig groups handle the index incrementations
		// by themselves
		// TODO: create the records iterator and pass it
		if (flag == ReadContigGroupStart)
		{
			DecompressConsensus(reads_, recIdx);
		}
		else if (flag == ReadTreeGroupStart)
		{
			DecompressTree(reads_[recIdx - 1],
						   reads_,
						   recIdx);
		}
		else if (flag == ReadIdentical)
		{
			ASSERT(lzContextStack.size() > 0);
			LzContext& lzCtx = lzContextStack.top();

			const FastqRecord* match = *(lzCtx.history.rbegin());

			DecompressRecord(reads_[recIdx], (ReadFlags)flag, match);
			recIdx++;
		}
		else
		{
			ASSERT(flag < ReadTreeGroupStart);
			DecompressRecord(reads_[recIdx], (ReadFlags)flag);
			recIdx++;
		}
	}
}


void LzDecompressorSE::DecompressLzMatch(FastqRecord& rec_, uint32 flag_)
{
	ASSERT(lzContextStack.size() > 0);
	LzContext& lzCtx = lzContextStack.top();

	rec_.seqLen = blockDesc.header.recMinLen;

	ASSERT(flag_ == ReadFullEncode || flag_ == ReadShiftOnly || flag_ == ReadFullExpensive);


	ASSERT(rec_.seqLen > 0);
	ASSERT(rec_.seqLen < 255);

	int32 recLen = rec_.seqLen;
	char* recStr = rec_.seq;

	int32 shift = (int32)(mainCtx.readers[FastqWorkBuffersSE::ShiftBuffer]->GetByte()) - ShiftOffset;

	ASSERT(ABS(shift) < recLen - params.minimizer.signatureLen);

	int32 prevId = mainCtx.dna.lzRle0Coder->GetSym();
	const FastqRecord* bestRec = (const FastqRecord*)*(lzCtx.history.rbegin() + prevId);

	const char* bestSeq = bestRec->seq;
	int32 bestLen = bestRec->seqLen;
	int32 bestPos = bestRec->minimPos;

	int32 nextMinPos = 0;

	if (shift >= 0)
	{
		bestSeq += shift;
		bestLen -= shift;
		bestPos -= shift;
		nextMinPos = bestPos;
	}
	else
	{
		for (int32 i = 0; i < -shift; ++i)
		{
			int32 c = mainCtx.dna.lettersXRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int)'N']);
			ASSERT(c < 5);
			c = idxToDna[c];

			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			rec_.seq[i] = c;
		}

		recStr += -shift;
		recLen -= -shift;
		nextMinPos = bestPos + (-shift);
	}

	// the '<=' was changed from '<' due to compression of rev-compl whole tr
	ASSERT(nextMinPos < rec_.seqLen);

	rec_.minimPos = nextMinPos;

	int32 minLen = MIN(bestLen, recLen);
	ASSERT(minLen > 0);

	if (flag_ == ReadFullEncode)
	{
		// TODO: optimize this routine
		for (int32 i = 0; i < minLen; ++i)
		{
			if (i == bestPos)
			{
				std::copy(lzCtx.currentMinimizerBuf.data(),
						  lzCtx.currentMinimizerBuf.data() + params.minimizer.signatureLen,
						  recStr + i);
				i += params.minimizer.signatureLen - 1;
				continue;
			}

			if (mainCtx.dna.matchRleCoder->GetSym())
			{
				recStr[i] = bestSeq[i];
			}
			else
			{
				int32 c = mainCtx.dna.lettersXRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int32)bestSeq[i]]);
				ASSERT(c < 5);
				c = idxToDna[c];

				ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
				recStr[i] = c;
			}
		}
	}
	else if (flag_ == ReadFullExpensive)
	{
		for (int32 i = 0; i < minLen; ++i)
		{
			if (i == bestPos)
			{
				std::copy(lzCtx.currentMinimizerBuf.data(),
						  lzCtx.currentMinimizerBuf.data() + params.minimizer.signatureLen,
						  recStr + i);
				i += params.minimizer.signatureLen - 1;
				continue;
			}

			if (mainCtx.dna.matchRcCoder->coder.DecodeSymbol(mainCtx.dna.matchRcCoder->rc) != 0)
			{
				recStr[i] = bestSeq[i];
			}
			else
			{
				int32 c = mainCtx.dna.lettersXRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int32)bestSeq[i]]);
				ASSERT(c < 5);
				c = idxToDna[c];

				ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
				recStr[i] = c;
			}
		}
	}
	else					// all symbols match
	{
		// TODO: use std::copy
		for (int32 i = 0; i < minLen; ++i)
			recStr[i] = bestSeq[i];
	}

	for (int32 i = minLen; i < recLen; ++i)
	{
		int32 c = mainCtx.dna.lettersXRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersXRcCoder->rc, dnaToIdx[(int)'N']);
		ASSERT(c < 5);
		c = idxToDna[c];

		ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
		recStr[i] = c;
	}

	bool isRev = mainCtx.dna.revRcCoder->coder.DecodeSymbol(mainCtx.dna.revRcCoder->rc) != 0;
	rec_.SetReadReverse(isRev);

	// update LZ context
	lzCtx.history.push_back((const FastqRecord*)&rec_);
}


void LzDecompressorSE::DecompressHardRead(FastqRecord& rec_)
{
	ASSERT(lzContextStack.size() > 0);
	LzContext& lzCtx = lzContextStack.top();

	rec_.seqLen = blockDesc.header.recMinLen;
	rec_.SetReadReverse(mainCtx.dna.revRcCoder->coder.DecodeSymbol(mainCtx.dna.revRcCoder->rc) != 0);

	ASSERT(rec_.seqLen > 0);
	ASSERT(rec_.seqLen < 255);

	rec_.minimPos = 0;
	for (int32 i = 0; i < rec_.seqLen; ++i)
	{
#if (ENC_HR_AC)
		uint32 x = mainCtx.dna.hardCoder->coder.DecodeSymbol(mainCtx.dna.hardCoder->rc);
		ASSERT(x <= 5);

		if (x < 5)
		{
			char c = (char)idxToDna[x];
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			rec_.seq[i] = c;
		}
		else
		{
			rec_.minimPos = i;
			std::copy(lzCtx.currentMinimizerBuf.data(),
					  lzCtx.currentMinimizerBuf.data() + params.minimizer.signatureLen,
					  rec_.seq + i);
			i += params.minimizer.signatureLen - 1;
		}
#else
		int32 c = mainCtx.readers[FastqWorkBuffersSE::HardReadsBuffer]->GetByte();				// TODO: read and copy N bytes
		if (c != MinimizerPositionSymbol)
		{
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			rec_.seq[i] = (char)c;
		}
		else
		{
			rec_.minimPos = i;
			std::copy(lzCtx.currentMinimizerBuf.data(),
					  lzCtx.currentMinimizerBuf.data() + params.minimizer.signatureLen,
					  rec_.seq + i);
			i += params.minimizer.signatureLen - 1;
		}
#endif
	}

	// update lz context
	//
	lzCtx.history.clear();
	lzCtx.history.push_back((const FastqRecord*)&rec_);
}


void LzDecompressorSE::DecompressExactRead(FastqRecord& rec_, const FastqRecord& aux_)
{
	// the exact match is identical as the previous read - only a different
	// direction is possible
	//
	rec_.SetReadReverse(mainCtx.dna.revRcCoder->coder.DecodeSymbol(mainCtx.dna.revRcCoder->rc) != 0);

	std::copy(aux_.seq, aux_.seq + aux_.seqLen, rec_.seq);
	rec_.seqLen = aux_.seqLen;
	rec_.minimPos = aux_.minimPos;
}


// TODO: this one should be the same as the upper one
void LzDecompressorSE::DecompressConsensusExactMatch(FastqRecord& rec_, const FastqRecord& prev_)
{
	rec_.SetReadReverse(mainCtx.dna.revRcCoder->coder.DecodeSymbol(mainCtx.dna.revRcCoder->rc) != 0);

	std::copy(prev_.seq, prev_.seq + prev_.seqLen, rec_.seq);
	rec_.seqLen = prev_.seqLen;
	rec_.minimPos = prev_.minimPos;
}


void LzDecompressorSE::DecompressConsensus(std::vector<FastqRecord>& reads_, uint64& startIdx_)
{
	// initialize lz context and main consensus encoding read
	//
	ASSERT(lzContextStack.size()  > 0);
	LzContext& lzCtx = lzContextStack.top();

	const FastqRecord* mainRec = *(lzCtx.history.rbegin());


	// initialize consensus context
	//
	ASSERT(blockDesc.isLenConst);
	const uint32 seqLen = blockDesc.header.recMinLen;			// TODO: MAX LEN

	consDecoders.push(ConsensusDecoder(seqLen));
	ConsensusDecoder& consDecoder = consDecoders.top();

	consDecoder.Reset(seqLen);
	consDecoder.type = ConsensusDecoder::ConsensusGroup;
	consDecoder.lastMinimPos = mainRec->minimPos;

	std::copy(lzCtx.currentMinimizerBuf.data(),
			  lzCtx.currentMinimizerBuf.data() + params.minimizer.signatureLen,
			  consDecoder.curSignatureString.data());

	if (seqLen < 128)
	{
		const int32 r1Rescale = ShiftOffset - consDecoder.seqLen / 2;
		const int32 r2Rescale = ShiftOffset - consDecoder.seqLen * 3 / 2;
		consDecoder.definition.range.first = (uint32)((int32)mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte() - r1Rescale);
		consDecoder.definition.range.second = (uint32)((int32)mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte() - r2Rescale);
	}
	else
	{
		consDecoder.definition.range.first = (uint32)(mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte());

		//int32 rangeDiff = mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte();

		// delta = range.second - range.first - seqLen + rescale
		int32 delta = mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte();
		int32 rescale = params.minimizer.signatureLen + ReadsContigBuilderParams::Default::BeginCut + ReadsContigBuilderParams::Default::EndCut;

		int32 rangeDiff = delta - rescale + seqLen;

		consDecoder.definition.range.second = consDecoder.definition.range.first + rangeDiff;

		ASSERT(consDecoder.definition.range.second < seqLen * 2);
	}


	// calculate the range of lz-match
	//
	std::pair<uint32, uint32> refRange;
	refRange.first = seqLen - mainRec->minimPos;
	refRange.second = refRange.first + seqLen;
	const uint32 rangeSize = consDecoder.definition.range.second - consDecoder.definition.range.first;

	for (uint32 i = consDecoder.definition.range.first; i < consDecoder.definition.range.first + rangeSize; ++i)
	{
		// skip the signature position range
		//
		if (i == seqLen)
		{
			i += params.minimizer.signatureLen - 1;
			continue;
		}

		consDecoder.definition.variantPositions[i] = !mainCtx.dna.consMatchCoder->GetSym();

		// read cons symbols both when there's a mismatch -- for letters context
		if (i < refRange.first + ReadsContigBuilderParams::Default::BeginCut
				|| i >= refRange.second - ReadsContigBuilderParams::Default::EndCut
				|| consDecoder.definition.variantPositions[i])
		{
			int32 c = mainCtx.dna.lettersCRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersCRcCoder->rc, dnaToIdx[(int)'N']);
			ASSERT(c < 5);
			consDecoder.definition.sequence[i] = idxToDna[c];
		}
		else
		{
			//SOFT_ASSERT(mainRec->seq[i - refRange.first] != 'N');
			consDecoder.definition.sequence[i] = mainRec->seq[i - refRange.first];
		}
	}
	consDecoder.readsCount = 0;


	// decompress first record of consensus without reading next flag
	//
	bool firstRead = true;
	uint32 decompTrees = 0;
	while (mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->PeekByte() != ReadGroupEnd)
	{
		ASSERT(startIdx_ < reads_.size());

		uint32 flag;
		if (!firstRead)
			flag = mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->GetByte();
		else
			flag = ReadContigGroupNext;

		ASSERT(flag == ReadContigGroupNext || flag == ReadTreeGroupStart || flag == ReadIdentical);

		if (flag == ReadIdentical)
		{
			FastqRecord& rec = reads_[startIdx_];
			DecompressRecord(rec, ReadIdentical, &reads_[startIdx_ - 1]);
			startIdx_++;
		}
		else if (flag == ReadTreeGroupStart)
		{
			DecompressTree(reads_[startIdx_ - 1], reads_, startIdx_);
			decompTrees++;
		}
		else
		{
			FastqRecord& rec = reads_[startIdx_++];
			DecompressRecord(rec, ReadContigGroupNext, NULL, firstRead);
			firstRead = false;

			// update lz buffer -- here we canot reuse the initial lzCtx local variable
			// as the lzContextStack could have been resized in nested method
			// invalidating all previous refenrences to its elements
			LzContext& lzCtx = lzContextStack.top();
			lzCtx.history.push_back((const FastqRecord*)&rec);
		}
	}

	// assert that we fininshed decoding of the contig
	//
	uint32 flag = mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->GetByte();

	ASSERT(flag == ReadGroupEnd);
	ASSERT(mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->PeekByte() != ReadContigGroupStart
			&& mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->PeekByte() != ReadTreeGroupStart);


	// clear the contig context
	//
	consDecoders.pop();
}


#include <stdio.h>


void LzDecompressorSE::DecompressConsensusRead(FastqRecord& rec_, bool useTreeShiftBuffer_)
{
	// decompress rev-compl
	//
	rec_.SetReadReverse(mainCtx.dna.revRcCoder->coder.DecodeSymbol(mainCtx.dna.revRcCoder->rc) != 0);

	// WARN: assume constant length
	rec_.seqLen = blockDesc.header.recMinLen;


	ConsensusDecoder& consDecoder = consDecoders.top();

#if 0
	// decode signature pos -- using delta encoding
	//
	if (!useTreeShiftBuffer_)
	{
		// TODO: handle case when readLen > 255
		//

		int32 dpos = (int32)mainCtx.readers[FastqWorkBuffersSE::ConsensusShiftBuffer]->GetByte() - ShiftOffset;
		ASSERT((int32)consDecoder.lastMinimPos + dpos < (int32)consDecoder.seqLen);
		ASSERT((int32)consDecoder.lastMinimPos + dpos >= 0);

		rec_.minimPos = (uint16)(consDecoder.lastMinimPos + dpos);
		ASSERT(rec_.minimPos < rec_.seqLen);
	}
	else
	{
		if (rec_.seqLen * 2 >= 256)
		{
			if (rec_.seqLen < 256)
			{
				rec_.minimPos = (int32)mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte();
				ASSERT(rec_.minimPos < rec_.seqLen);
			}
			else
			{
				/* TODO: special case for variable length reads
				int32 dpos = (int32)mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->Get2Bytes();
				dpos -= (int32)rec_.seqLen;

				ASSERT((int32)consDecoder.lastMinimPos + dpos < (int32)consDecoder.seqLen);

				rec_.minimPos = (uint16)(consDecoder.lastMinimPos + dpos);
				ASSERT(rec_.minimPos < rec_.seqLen);
				*/
				rec_.minimPos = (int32)mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->Get2Bytes();
				ASSERT(rec_.minimPos < rec_.seqLen);
			}
		}
		else
		{
			int32 dpos = (int32)mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte() - ShiftOffset;

			ASSERT(consDecoder.lastMinimPos + dpos >= 0);
			ASSERT((int32)consDecoder.lastMinimPos + dpos < (int32)consDecoder.seqLen);

			rec_.minimPos = (uint16)(consDecoder.lastMinimPos + dpos);
			ASSERT(rec_.minimPos < rec_.seqLen);
		}
	}
#else

	BitMemoryReader* reader = NULL;
	if (!useTreeShiftBuffer_)
	{
		reader = mainCtx.readers[FastqWorkBuffersSE::ConsensusShiftBuffer];
	}
	else
	{
		reader = mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer];
	}

	if (rec_.seqLen * 2 >= 256)
	{
		ASSERT(rec_.seqLen < 256);

		rec_.minimPos = reader->GetByte();
		ASSERT(rec_.minimPos < rec_.seqLen);
	}
	else
	{
		int32 dpos = (int32)reader->GetByte() - ShiftOffset;

		ASSERT(consDecoder.lastMinimPos + dpos >= 0);
		ASSERT((int32)consDecoder.lastMinimPos + dpos < (int32)consDecoder.seqLen);

		rec_.minimPos = (uint16)(consDecoder.lastMinimPos + dpos);
		ASSERT(rec_.minimPos < rec_.seqLen);

	}


#endif

	consDecoder.lastMinimPos = rec_.minimPos;


	// calculate the read range sequence inside consensus
	//
	ASSERT(consDecoder.seqLen >= (uint32)rec_.minimPos || !printf("%d %d\n", consDecoder.seqLen, rec_.minimPos));
	const uint32 consStart = (uint32)consDecoder.seqLen - (uint32)rec_.minimPos;

	// TODO: simplify
	ASSERT(consStart + ReadsContigBuilderParams::Default::BeginCut >= consDecoder.definition.range.first
		   || (rec_.minimPos <= ReadsContigBuilderParams::Default::BeginCut
			   && consStart + rec_.minimPos + params.minimizer.signatureLen >= consDecoder.definition.range.first)
		   || printf("WARN: (%d) %d %d %d\n", blockDesc.header.minimizerId, consStart, consDecoder.definition.range.first, rec_.minimPos) == 0);
	ASSERT(consStart + consDecoder.seqLen - ReadsContigBuilderParams::Default::EndCut <= consDecoder.definition.range.second
		   || ((int32)rec_.seqLen - rec_.minimPos - params.minimizer.signatureLen <= (int32)ReadsContigBuilderParams::Default::EndCut
			   && consStart + rec_.minimPos - ReadsContigBuilderParams::Default::EndCut <= consDecoder.definition.range.second)
		   || printf("WARN: min pos: %d ; cons start: %d ; cons end: %d \n",
									  rec_.minimPos,
									  consDecoder.definition.range.first,
									  consDecoder.definition.range.second) == 0);


	// read the begin cut of the read
	//
	uint32 decodeIterator = 0;
	while (decodeIterator < ReadsContigBuilderParams::Default::BeginCut)
	{
		if (decodeIterator == rec_.minimPos)
		{
			std::copy(consDecoder.curSignatureString.data(),
					  consDecoder.curSignatureString.data() + params.minimizer.signatureLen,
					  rec_.seq + decodeIterator);

			ASSERT(decodeIterator + params.minimizer.signatureLen > ReadsContigBuilderParams::Default::BeginCut);
			decodeIterator += params.minimizer.signatureLen;
			continue;
		}

		ASSERT(consStart + decodeIterator < consDecoder.definition.sequence.size()
			   || !printf("%d %d %d\n", consStart, decodeIterator, (int32)consDecoder.definition.sequence.size()));
		char ctx = consDecoder.definition.sequence[consStart + decodeIterator];
		ASSERT(ctx == 'A' || ctx == 'C' || ctx == 'G' || ctx == 'T' || ctx == 'N');
		int32 c = mainCtx.dna.lettersCRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersCRcCoder->rc,
													  dnaToIdx[(int)consDecoder.definition.sequence[consStart + decodeIterator]]);
		ASSERT(c < 5);
		rec_.seq[decodeIterator] = idxToDna[c];

		decodeIterator++;
	}


	// read the main read part
	//
	while (decodeIterator < consDecoder.seqLen - ReadsContigBuilderParams::Default::EndCut)
	{
		// skip reading of the signature position range
		//
		if (decodeIterator == rec_.minimPos)
		{
			std::copy(consDecoder.curSignatureString.data(),
					  consDecoder.curSignatureString.data() + params.minimizer.signatureLen,
					  rec_.seq + decodeIterator);

			decodeIterator += params.minimizer.signatureLen;
			continue;
		}

		ASSERT(consStart + decodeIterator >= consDecoder.definition.range.first);

		//ASSERT(consDecoder.sequence[consStart + j] != '.');
		if (consDecoder.definition.variantPositions[consStart + decodeIterator])
		{
			int32 idx = mainCtx.dna.lettersCRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersCRcCoder->rc,
															dnaToIdx[(int)consDecoder.definition.sequence[consStart + decodeIterator]]);
			ASSERT(idx < 5);
			rec_.seq[decodeIterator] = idxToDna[idx];
		}
		else
		{
			//ASSERT(consDecoder.definition.sequence[consStart + decodeIterator] != 'N');
			rec_.seq[decodeIterator] = consDecoder.definition.sequence[consStart + decodeIterator];
		}

		decodeIterator++;
	}


	// read the end cut
	//
	while (decodeIterator < consDecoder.seqLen)
	{
		// end zone + end cut
		// this condition will be false in case of miniPos == seq_len - signatureLen - [0, end_cut]
		//ASSERT(consStart + decodeIterator <= consDecoder.definition.range.second + ReadsContigBuilderParams::Default::EndCut);

		int32 c = mainCtx.dna.lettersCRcCoder->coder.DecodeSymbol(mainCtx.dna.lettersCRcCoder->rc,
													  dnaToIdx[(int)consDecoder.definition.sequence[consStart + decodeIterator]]);
		ASSERT(c < 5);
		rec_.seq[decodeIterator] = idxToDna[c];

		decodeIterator++;
	}

	consDecoder.readsCount++;
}


void LzDecompressorSE::DecompressTree(const FastqRecord& rootRec_, std::vector<FastqRecord>& reads_, uint64& startIdx_)
{
	// create a local root version of main encoding read
	//
	FastqRecord mainLzRec = rootRec_;
	uint32 altMinimPos = 0;

	if (rootRec_.seqLen * 2 >= 256)
	{
		ASSERT(rootRec_.seqLen < 256);

		altMinimPos = mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte();
	}
	else
	{
		// shift = mainTree.minimizerPos - rootNode.minimizerPos
		int32 shift = (int32)mainCtx.readers[FastqWorkBuffersSE::TreeShiftBuffer]->GetByte() - ShiftOffset;

		altMinimPos = shift + (int32)mainLzRec.minimPos;
	}

	ASSERT(altMinimPos <= (uint32)mainLzRec.seqLen - (uint32)params.minimizer.signatureLen);

	mainLzRec.minimPos = altMinimPos;


	// create new lz context
	//
	lzContextStack.push(LzContext());
	LzContext& curLzCtx = lzContextStack.top();
	curLzCtx.history.push_back((const FastqRecord*)&mainLzRec);

	std::copy(mainLzRec.seq + altMinimPos,
			  mainLzRec.seq + altMinimPos + params.minimizer.signatureLen,
			  curLzCtx.currentMinimizerBuf.data());


	// decode reads
	//
	while (mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->PeekByte() != ReadGroupEnd)
	{
		uint32 flag = mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->GetByte();
		//ASSERT(flag != ReadDifficult);

		if (flag == ReadContigGroupStart)
		{
			DecompressConsensus(reads_, startIdx_);
		}
		else if (flag == ReadTreeGroupStart)
		{
			DecompressTree(reads_[startIdx_ - 1],
						   reads_,
						   startIdx_);
		}
		else if (flag == ReadIdentical)
		{
			ASSERT(lzContextStack.size() > 0);
			LzContext& lzCtx = lzContextStack.top();

			const FastqRecord* match = *(lzCtx.history.rbegin());

			DecompressRecord(reads_[startIdx_], (ReadFlags)flag, match);
			startIdx_++;
		}
		else
		{
			ASSERT(flag < ReadContigGroupStart);
			DecompressRecord(reads_[startIdx_], (ReadFlags)flag);
			startIdx_++;
		}
	}


	// assert the end of grop decoding
	//
	uint32 flag = mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->GetByte();
	ASSERT(flag == ReadGroupEnd);


	// handle the multi-tree case !!!
	//
	if (mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->PeekByte() == ReadTreeGroupStart)
	{
		mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->GetByte();

		DecompressTree(rootRec_, reads_, startIdx_);
	}


	// clean decoders
	//
	lzContextStack.pop();
}


void LzDecompressorSE::DecompressRecord(FastqRecord& rec_, ReadFlags flag_,
										const FastqRecord* auxRec_, bool aux_)
{
	// decompress header
	//


	// decompress dna data
	//
	switch (flag_)
	{
	case ReadDifficult:
		DecompressHardRead(rec_);
		break;

	case ReadIdentical:
	{
		ASSERT(auxRec_ != NULL);
		DecompressExactRead(rec_, *auxRec_);
		break;
	}

	case ReadShiftOnly:
	case ReadFullEncode:
	case ReadFullExpensive:
		DecompressLzMatch(rec_, flag_);
		break;

	case ReadContigGroupNext:
		DecompressConsensusRead(rec_, aux_);
		break;

	default: ASSERT(0);
	}


	// decompress quality data
	//
	DecompressReadQuality(rec_);
}

void LzDecompressorSE::DecompressMetaAndLinkRecords(std::vector<FastqRecord> &reads_, FastqChunk &chunk_)
{
	// TODO: now only supports const-length records
	//
	ASSERT(blockDesc.header.recMinLen == blockDesc.header.recMaxLen);

	// TODO: decompress read lengths here
	//

	char* dnaBufferPtr = (char*)chunk_.data.Pointer();
	uint64 pos = 0;
	for (uint64 i = 0; i < reads_.size(); ++i)
	{
		FastqRecord& r = reads_[i];
		r.seqLen = blockDesc.header.recMaxLen;

		// decompress header
		//
		if (params.archType.readsHaveHeaders)
		{
			r.head = dnaBufferPtr + pos;
			DecompressReadId(r, mainCtx.id);
			pos += r.headLen;
		}

		// decompress read length
		//
		r.seq = dnaBufferPtr + pos;
		r.qua = r.seq + r.seqLen;

		pos += r.seqLen * 2;

	}

	chunk_.size = pos;
}

void LzDecompressorSE::DecompressReadQuality(FastqRecord &rec_)
{
	IQualityStoreBase::DecompressReadQuality(rec_, mainCtx.qua, params);
}

void LzDecompressorSE::DecompressQuality()
{
	// HINT: here can be performed extra quality information compressed
	//
	// ...
	//
}



// RawCompressorSE / RawDecompressorSE
//
//
RawCompressorSE::RawCompressorSE(const CompressorParams& params_,
								 const QualityCompressionData& globalQuaData_,
								 const FastqRawBlockStats::HeaderStats& headData_,
								 const CompressorAuxParams& auxParams_)
	:	IDnaRawStoreBase(params_, globalQuaData_, headData_, auxParams_)
{
	ppmdCoder = new PpmdEncoder();
	ppmdCoder->StartCompress(DefaultPpmdOrder, DefaultPpmdMemorySizeMb);
}


RawCompressorSE::~RawCompressorSE()
{
	ppmdCoder->FinishCompress();
	delete ppmdCoder;
}


void RawCompressorSE::Compress(const std::vector<FastqRecord>& reads_,
									 PackContext& packCtx_,
									 uint32 minimizerId_,
									 uint64 rawDnaStreamSize_,
									 FastqCompressedBin& fastqWorkBin_,
									 CompressedFastqBlock &compBin_)
{
	ASSERT(reads_.size() != 0);
	ASSERT(packCtx_.stats.maxSeqLen > 0);
	ASSERT(packCtx_.stats.maxSeqLen >= packCtx_.stats.minSeqLen);


	// initialize block description header
	//
	blockDesc.Reset();
	blockDesc.header.minimizerId = minimizerId_;
	blockDesc.header.recordsCount = reads_.size();
	blockDesc.header.recMinLen = packCtx_.stats.minSeqLen;
	blockDesc.header.recMaxLen = packCtx_.stats.maxSeqLen;
	blockDesc.header.rawDnaStreamSize = rawDnaStreamSize_;
	blockDesc.isLenConst = (packCtx_.stats.minSeqLen == packCtx_.stats.maxSeqLen);

	// prepare the work and output buffer buffer
	//
	const uint32 bitsPerLen = blockDesc.isLenConst ? 0 : bit_length(blockDesc.header.recMaxLen - blockDesc.header.recMinLen);
	const uint64 dnaPreallocSize = rawDnaStreamSize_ * 2;			// account here for both DNA + quality

	uint64 idPreallocSize = 0;
	if (params.archType.readsHaveHeaders)
	{
		// estimate prealloc size
		const auto& r = reads_.front();
		idPreallocSize = (r.headLen * 1.5) * blockDesc.header.recordsCount;
	}

	const uint64 compPreallocSize = dnaPreallocSize + idPreallocSize + (RawBlockHeader::Size)
									+ (((bitsPerLen > 0) ? blockDesc.header.recordsCount * bitsPerLen : 0));


	// preapre buffers and readers
	//
	DataChunk& dnaWorkBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::HardReadsBuffer];
	DataChunk& quaWorkBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::QualityBuffer];

	if (dnaWorkBuffer.data.Size() < dnaPreallocSize)
		dnaWorkBuffer.data.Extend(dnaPreallocSize);
	dnaWorkBuffer.size = 0;

	if (quaWorkBuffer.data.Size() < dnaPreallocSize)
		quaWorkBuffer.data.Extend(dnaPreallocSize);
	quaWorkBuffer.size = 0;

	if (compBin_.dataBuffer.data.Size() < compPreallocSize)
		compBin_.dataBuffer.data.Extend(compPreallocSize);
	compBin_.dataBuffer.size = 0;


	// start encoding
	//
	StartEncoding(fastqWorkBin_.buffers);


	// sort records
	//
	TFastqComparator<const MatchNode&> comparator;
	std::sort(packCtx_.graph->nodes.begin(), packCtx_.graph->nodes.end(), comparator);


	ASSERT(blockDesc.isLenConst);


	// compress records
	//
	for (const FastqRecord& rec : reads_)
	{
		if (auxParams.dry_run)
		{
			DryRunCompress(rec, dryFastqWriter);
		}
		else
		{
			// compress meta data (e.g. read length) and read id
			//
			if (params.archType.readsHaveHeaders)
			{
				CompressReadId(rec, mainCtx.id);
				blockDesc.header.rawIdStreamSize += rec.headLen;
			}

			CompressReadSequence(rec);

			CompressReadQuality(rec);
		}
	}

	// end encoding
	//
	EndEncoding(fastqWorkBin_.buffers);



	// prepare output writer and reserve header space
	//
	BitMemoryWriter blockWriter(compBin_.dataBuffer.data);
	blockWriter.FillBytes(0, RawBlockHeader::Size);


	// compress header buffers (of present)
	//
	// -- now just copy, will be done with more advanced codecs --

	uint64 headOffset = blockWriter.Position();
	if (params.archType.readsHaveHeaders)
	{
		DataChunk& tokenBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::ReadIdTokenBuffer];
		DataChunk& valueBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::ReadIdValueBuffer];

		std::copy(tokenBuffer.data.Pointer(),
				  tokenBuffer.data.Pointer() + tokenBuffer.size,
				  compBin_.dataBuffer.data.Pointer() + headOffset);

		std::copy(valueBuffer.data.Pointer(),
				  valueBuffer.data.Pointer() + valueBuffer.size,
				  compBin_.dataBuffer.data.Pointer() + headOffset + tokenBuffer.size);

		blockDesc.header.idTokenCompSize = tokenBuffer.size;
		blockDesc.header.idValueCompSize = valueBuffer.size;
	}
	else
	{
		blockDesc.header.idTokenCompSize = 0;
		blockDesc.header.idValueCompSize = 0;
	}



	// compress dna
	//
	const uint64 dnaDataOffset = headOffset + blockDesc.header.idTokenCompSize + blockDesc.header.idValueCompSize;
	uint64 compDnaSize = compBin_.dataBuffer.data.Size() - dnaDataOffset;

	if (!auxParams.dry_run)
	{
		CompressDna(dnaWorkBuffer, dnaWorkBuffer.size, compBin_.dataBuffer, compDnaSize, dnaDataOffset);
		ASSERT(compDnaSize > 0);
	}

	blockDesc.header.dnaCompSize = compDnaSize;


	// compress quality
	//
	const uint64 quaDataOffset = dnaDataOffset + compDnaSize;
	uint64 compQuaSize = compBin_.dataBuffer.data.Size() - quaDataOffset;

	if (!auxParams.dry_run)
	{
		CompressQuality(quaWorkBuffer, quaWorkBuffer.size, compBin_.dataBuffer, compQuaSize, quaDataOffset);
		ASSERT(compQuaSize > 0);
	}

	blockDesc.header.quaCompSize = compQuaSize;


	// store footer and header
	//
	uint64 endPos = quaDataOffset + compQuaSize;

	blockWriter.SetPosition(endPos);
	StoreRawFooter(blockDesc.footer, blockWriter);

	blockDesc.header.footerOffset = endPos;
	blockDesc.header.footerSize = blockWriter.Position() - blockDesc.header.footerOffset;


	StoreHeader(blockDesc.header, compBin_.dataBuffer);

	compBin_.dataBuffer.size = dnaDataOffset
			+ blockDesc.header.dnaCompSize + blockDesc.header.quaCompSize
			+ blockDesc.header.footerSize;
	ASSERT(compBin_.dataBuffer.size == blockDesc.header.footerOffset + blockDesc.header.footerSize);


	// save for debug purposes
	//
	compBin_.stats.counts["NDnaCompSize"] = blockDesc.header.dnaCompSize;
	compBin_.stats.counts["NQuaCompSize"] = blockDesc.header.quaCompSize;
	compBin_.stats.counts["NDnaRawSize"] = blockDesc.header.rawDnaStreamSize;
	compBin_.stats.counts["NReadIdTokenCompSize"] = blockDesc.header.idTokenCompSize;
	compBin_.stats.counts["NReadIdValueCompSize"] = blockDesc.header.idValueCompSize;
	compBin_.stats.counts["NReadIdRawSize"] = blockDesc.header.rawIdStreamSize;

	compBin_.signatureId = minimizerId_;
}


void RawCompressorSE::DryRunCompress(const FastqRecord &record_, BitMemoryWriter* dryFastqWriter_)
{
	// copied from LzCompressor -- TODO: merge into one method to de-duplicate the code
	//
	ASSERT(auxParams.f_uncompressed != NULL);

	if (auxParams.output_fastq)
	{
		// output reads header
		//
		if (params.archType.readsHaveHeaders)
		{
			dryFastqWriter_->PutBytes((byte*)record_.head, record_.headLen);
		}
		else // dummy header -- this should not be the case
		{
			//const std::string libName = "SRX000000.X";
			//dryFastqWriter_->PutBytes((byte*)libName.c_str(), libName.size());

			char buf[64];
			int len = 0;
			buf[len++] = '@';
			params.minimizer.GenerateMinimizer(params.minimizer.TotalMinimizersCount(), buf + len);
			len += params.minimizer.signatureLen;
			buf[len++] = '.';
			len += to_string(buf + len, currentRecordIdx);

			dryFastqWriter_->PutBytes((byte*)buf, len);

		}
		dryFastqWriter_->PutByte('\n');


		// output dna
		//
		if (record_.IsReadReverse())
		{
			const char* codes = FastqRecord::GetRCCodes();
			for (int32 i = record_.seqLen - 1; i >= 0; i--)
				dryFastqWriter_->PutByte(codes[record_.seq[i] - 64]);
		}
		else
		{
			dryFastqWriter_->PutBytes((byte*)record_.seq, record_.seqLen);
		}
		dryFastqWriter_->PutByte('\n');

		// otput control line
		dryFastqWriter_->PutByte('+');
		dryFastqWriter_->PutByte('\n');

		currentRecordIdx++;
	}


	// output quality
	//
	{
		std::vector<byte> qscores;
		IQualityStoreBase::CompressReadQuality(record_, mainCtx.qua, params, true, &qscores);

		dryFastqWriter_->PutBytes(qscores.data(), qscores.size());
		dryFastqWriter_->PutByte('\n');
	}
}


void RawCompressorSE::StartEncoding(std::vector<DataChunk*>& buffers_)
{
	// create main readers
	//
	mainCtx.dnaWriter = new BitMemoryWriter(buffers_[FastqWorkBuffersSE::HardReadsBuffer]->data);
	mainCtx.quaWriter = new BitMemoryWriter(buffers_[FastqWorkBuffersSE::QualityBuffer]->data);


	// link readers or create encoder decorators
	//
	mainCtx.dna.rawCoder = mainCtx.dnaWriter;

	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		mainCtx.qua.rawCoder = mainCtx.quaWriter;
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		mainCtx.qua.binaryCoder = new QualityEncoders::BinaryEncoder(*mainCtx.quaWriter);
		mainCtx.qua.binaryCoder->Start();
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		mainCtx.qua.illu8Coder = new QualityEncoders::Illu8Encoder(*mainCtx.quaWriter);
		mainCtx.qua.illu8Coder->Start();
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		mainCtx.qua.qvzCoder = new /*QualityEncoders::*/QVZEncoder(mainCtx.quaWriter, globalQuaData.codebook.qlist);
		mainCtx.qua.qvzCoder->Start();


		// also reset rng when starting encoding so that the random number
		// generator will have the same initial seed also when
		// compressing in multithreaded mode
		ResetWellRng();
		break;
	}
	}


	// create read id encoders
	//
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.idTokenWriter = new BitMemoryWriter(buffers_[FastqWorkBuffersSE::ReadIdTokenBuffer]->data);
		mainCtx.idValueWriter = new BitMemoryWriter(buffers_[FastqWorkBuffersSE::ReadIdValueBuffer]->data);

		mainCtx.id.tokenCoder = new FieldEncoders::TokenEncoder(*mainCtx.idTokenWriter);
		mainCtx.id.valueCoder = new FieldEncoders::ValueEncoder(*mainCtx.idValueWriter);

		mainCtx.id.tokenCoder->Start();
		mainCtx.id.valueCoder->Start();
	}


	// prepare writer for dry run
	//
	if (auxParams.dry_run)
	{
		dryFastqWriter->SetPosition(0);

		currentRecordIdx = 1;
	}
}


void RawCompressorSE::EndEncoding(std::vector<DataChunk*>& buffers_)
{
	// finish dna encoding
	//
	mainCtx.dnaWriter->FlushPartialWordBuffer();
	buffers_[FastqWorkBuffersSE::HardReadsBuffer]->size = mainCtx.dna.rawCoder->Position();

	mainCtx.dna.rawCoder = NULL;
	TFree(mainCtx.dnaWriter);


	// finish quality encoding
	//
	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		mainCtx.qua.rawCoder = NULL;		// this was just a pointer
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		mainCtx.qua.binaryCoder->End();
		TFree(mainCtx.qua.binaryCoder);
		break;
	}


	case QualityCompressionParams::MET_8BIN:
	{
		mainCtx.qua.illu8Coder->End();
		TFree(mainCtx.qua.illu8Coder);
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{	
		mainCtx.qua.qvzCoder->End();
		TFree(mainCtx.qua.qvzCoder);

		break;
	}
	}

	mainCtx.quaWriter->FlushPartialWordBuffer();
	buffers_[FastqWorkBuffersSE::QualityBuffer]->size = mainCtx.quaWriter->Position();

	TFree(mainCtx.quaWriter);


	// finish read id encoding
	//
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.id.tokenCoder->End();
		mainCtx.id.valueCoder->End();

		mainCtx.idTokenWriter->FlushPartialWordBuffer();
		mainCtx.idValueWriter->FlushPartialWordBuffer();

		buffers_[FastqWorkBuffersSE::ReadIdTokenBuffer]->size = mainCtx.idTokenWriter->Position();
		buffers_[FastqWorkBuffersSE::ReadIdValueBuffer]->size = mainCtx.idValueWriter->Position();

		TFree(mainCtx.id.tokenCoder);
		TFree(mainCtx.id.valueCoder);

		TFree(mainCtx.idTokenWriter);
		TFree(mainCtx.idValueWriter);
	}


	// finalize writer for dry run
	//
	if (auxParams.dry_run)
	{
		dryFastqWriter->Flush();

		// POSIX requires C stdio FILE* ops to be atomic -- no mutex here required
		//
		fwrite(dryFastqWriter->Pointer(), 1, dryFastqWriter->Position(), auxParams.f_uncompressed);
	}
}



void RawCompressorSE::CompressReadSequence(const FastqRecord &rec_)
{
	for (uint32 i = 0; i < rec_.seqLen; ++i)
		mainCtx.dna.rawCoder->PutByte(rec_.seq[i]);
}

void RawCompressorSE::CompressReadQuality(const FastqRecord &rec_)
{
	IQualityStoreBase::CompressReadQuality(rec_, mainCtx.qua, params);
}


void RawCompressorSE::CompressDna(const DataChunk &inChunk_, uint64 inSize_,
									 DataChunk &outChunk_, uint64 &outSize_, uint64 outOffset_)
{
	CompressBuffer(*ppmdCoder, inChunk_, inSize_, outChunk_, outSize_, outOffset_);
}

void RawCompressorSE::CompressQuality(const DataChunk &inChunk_, uint64 inSize_,
										 DataChunk &outChunk_, uint64 &outSize_, uint64 outOffset_)
{
	// HINT: here can be performed extra quality information compressed
	//
	// ...
	//


	// compress or copy buffers
	//
	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		CompressBuffer(*ppmdCoder, inChunk_, inSize_, outChunk_, outSize_, outOffset_);
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	case QualityCompressionParams::MET_8BIN:
	{
		// in case of using custom encoder, just copy the compressed binary data
		// to the output buffer
		std::copy(inChunk_.data.Pointer(), inChunk_.data.Pointer() + inSize_,
				  outChunk_.data.Pointer() + outOffset_);
		outSize_ = inSize_;
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{	
		// in case of using custom encoder, just copy the compressed binary data
		// to the output buffer
		std::copy(inChunk_.data.Pointer(), inChunk_.data.Pointer() + inSize_,
				  outChunk_.data.Pointer() + outOffset_);
		outSize_ = inSize_;

		break;
	}
	}
}



RawDecompressorSE::RawDecompressorSE(const CompressorParams& params_,
									 const QualityCompressionData& globalQuaData_,
									 const FastqRawBlockStats::HeaderStats& headData_)
	:	IDnaRawStoreBase(params_, globalQuaData_, headData_)
{
	ppmdCoder = new PpmdDecoder();
	ppmdCoder->StartDecompress(DefaultPpmdMemorySizeMb);
}


RawDecompressorSE::~RawDecompressorSE()
{
	ppmdCoder->FinishDecompress();
	delete ppmdCoder;
}


void RawDecompressorSE::Decompress(CompressedFastqBlock& compBin_,
										 std::vector<FastqRecord>& reads_,
										 FastqCompressedBin& fastqWorkBin_,
										 FastqChunk& dnaBuffer_)
{
	// read header and footer
	//
	blockDesc.Reset();
	ReadHeader(blockDesc.header, compBin_.dataBuffer);
	blockDesc.isLenConst = (blockDesc.header.recMinLen == blockDesc.header.recMaxLen);

	BitMemoryReader blockReader(compBin_.dataBuffer.data, compBin_.dataBuffer.size);
	blockReader.SetPosition(blockDesc.header.footerOffset);
	ReadRawFooter(blockDesc.footer, blockReader);
	blockReader.SetPosition(RawBlockHeader::Size);


	// prepare output dna bin
	//
	reads_.resize(blockDesc.header.recordsCount);

	uint64 bufSize = blockDesc.header.rawDnaStreamSize * 2;		// * 2 -- handle quality
	if (params.archType.readsHaveHeaders)
		bufSize += blockDesc.header.rawIdStreamSize;

	if (dnaBuffer_.data.Size() < bufSize)
		dnaBuffer_.data.Extend(bufSize);


	// preapre buffers
	//
	DataChunk& dnaWorkBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::HardReadsBuffer];
	DataChunk& quaWorkBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::QualityBuffer];

	if (dnaWorkBuffer.data.Size() < blockDesc.header.rawDnaStreamSize)
		dnaWorkBuffer.data.Extend(blockDesc.header.rawDnaStreamSize);

	if (quaWorkBuffer.data.Size() < blockDesc.header.rawDnaStreamSize)
		quaWorkBuffer.data.Extend(blockDesc.header.rawDnaStreamSize);


	// prepare output bin -- just raw reads
	//
	reads_.resize(blockDesc.header.recordsCount);


	// decompress meta data buffers and (if applicable) read id buffers
	//
	ASSERT(blockDesc.isLenConst);


	// decompress headers -- just copy the buffer contents
	//
	uint64 headOffset = blockReader.Position();
	if (params.archType.readsHaveHeaders)
	{
		DataChunk& tokenBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::ReadIdTokenBuffer];
		DataChunk& valueBuffer = *fastqWorkBin_.buffers[FastqWorkBuffersSE::ReadIdValueBuffer];

		if (tokenBuffer.data.Size() < blockDesc.header.idTokenCompSize)
			tokenBuffer.data.Extend(blockDesc.header.idTokenCompSize);

		if (valueBuffer.data.Size() < blockDesc.header.idValueCompSize)
			valueBuffer.data.Extend(blockDesc.header.idValueCompSize);

		std::copy(compBin_.dataBuffer.data.Pointer() + headOffset,
				  compBin_.dataBuffer.data.Pointer() + headOffset + blockDesc.header.idTokenCompSize,
				  tokenBuffer.data.Pointer());

		std::copy(compBin_.dataBuffer.data.Pointer() + headOffset + blockDesc.header.idTokenCompSize,
				  compBin_.dataBuffer.data.Pointer() + headOffset
					+ blockDesc.header.idTokenCompSize + blockDesc.header.idValueCompSize,
				  valueBuffer.data.Pointer());

		tokenBuffer.size = blockDesc.header.idTokenCompSize;
		valueBuffer.size = blockDesc.header.idValueCompSize;
	}



	// decompress dna
	//
	const uint64 dnaDataOffset = headOffset + blockDesc.header.idTokenCompSize + blockDesc.header.idValueCompSize;
	uint64 dnaOutSize = 0;
	DecompressDna(dnaWorkBuffer, dnaOutSize, compBin_.dataBuffer, blockDesc.header.dnaCompSize, dnaDataOffset);
	ASSERT(dnaOutSize == blockDesc.header.rawDnaStreamSize);
	dnaWorkBuffer.size = dnaOutSize;


	// decompress quality
	//
	const uint64 quaDataOffset = dnaDataOffset + blockDesc.header.dnaCompSize;
	uint64 quaOutSize = 0;
	DecompressQuality(quaWorkBuffer, quaOutSize, compBin_.dataBuffer, blockDesc.header.quaCompSize, quaDataOffset);
	quaWorkBuffer.size = quaOutSize;



	// start decoding
	//
	StartDecoding(fastqWorkBin_.buffers);


	// decode records
	//
	dnaBuffer_.size = blockDesc.header.rawDnaStreamSize * 2;
	if (params.archType.readsHaveHeaders)
		dnaBuffer_.size += blockDesc.header.rawIdStreamSize;

	DecompressRecords(reads_, dnaBuffer_);


	// end encoding
	//
	EndDecoding();
}


void RawDecompressorSE::StartDecoding(std::vector<DataChunk*>& buffers_)
{
	// start dna encoding
	//
	mainCtx.dnaReader = new BitMemoryReader(buffers_[FastqWorkBuffersSE::HardReadsBuffer]->data,
											buffers_[FastqWorkBuffersSE::HardReadsBuffer]->size);
	mainCtx.dna.rawCoder = mainCtx.dnaReader;		// just link


	// start quality encoding
	//
	mainCtx.quaReader = new BitMemoryReader(buffers_[FastqWorkBuffersSE::QualityBuffer]->data,
											buffers_[FastqWorkBuffersSE::QualityBuffer]->size);

	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		mainCtx.qua.rawCoder = mainCtx.quaReader;
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		mainCtx.qua.binaryCoder = new QualityDecoders::BinaryDecoder(*mainCtx.quaReader);
		mainCtx.qua.binaryCoder->Start();
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		mainCtx.qua.illu8Coder = new QualityDecoders::Illu8Decoder(*mainCtx.quaReader);
		mainCtx.qua.illu8Coder->Start();
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		mainCtx.qua.qvzCoder = new /*QualityDecoders::*/QVZDecoder(mainCtx.quaReader, globalQuaData.codebook.qlist);
		mainCtx.qua.qvzCoder->Start();


		// also reset rng when starting decoding so that the random number
		// generator will have the same initial seed also when
		// compressing in multithreaded mode
		ResetWellRng();
		break;
	}
	}


	// start read id encoding
	//
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.idTokenReader = new BitMemoryReader(buffers_[FastqWorkBuffersSE::ReadIdTokenBuffer]->data,
													buffers_[FastqWorkBuffersSE::ReadIdTokenBuffer]->size);
		mainCtx.idValueReader = new BitMemoryReader(buffers_[FastqWorkBuffersSE::ReadIdValueBuffer]->data,
													buffers_[FastqWorkBuffersSE::ReadIdValueBuffer]->size);

		mainCtx.id.tokenCoder = new FieldDecoders::TokenDecoder(*mainCtx.idTokenReader);
		mainCtx.id.valueCoder = new FieldDecoders::ValueDecoder(*mainCtx.idValueReader);

		mainCtx.id.tokenCoder->Start();
		mainCtx.id.valueCoder->Start();
	}
}


void RawDecompressorSE::EndDecoding()
{
	// finish dna encoding
	//
	mainCtx.dna.rawCoder = NULL;
	TFree(mainCtx.dnaReader);


	// finish quality decoding
	//
	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		mainCtx.qua.rawCoder = NULL;
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		mainCtx.qua.binaryCoder->End();
		TFree(mainCtx.qua.binaryCoder);
		break;
	}


	case QualityCompressionParams::MET_8BIN:
	{
		mainCtx.qua.illu8Coder->End();
		TFree(mainCtx.qua.illu8Coder);
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		mainCtx.qua.qvzCoder->End();
		TFree(mainCtx.qua.qvzCoder);

		break;
	}
	}

	TFree(mainCtx.quaReader);


	// fininsh headers decoders
	//
	if (params.archType.readsHaveHeaders)
	{
		mainCtx.id.tokenCoder->End();
		mainCtx.id.valueCoder->End();

		TFree(mainCtx.id.tokenCoder);
		TFree(mainCtx.id.valueCoder);

		TFree(mainCtx.idTokenReader);
		TFree(mainCtx.idValueReader);
	}
}


void RawDecompressorSE::DecompressRecords(std::vector<FastqRecord>& reads_, FastqChunk& buffer_)
{
	ASSERT(blockDesc.isLenConst);

	//const uint64 dnaBufferSize = blockDesc.header.rawDnaStreamSize;

	char* dnaBufferPtr = (char*)buffer_.data.Pointer();
	uint64 bufferPos = 0;

	for (FastqRecord& rec : reads_)
	{
		rec.Reset();


		// decompress meta data and read id
		//
		if (params.archType.readsHaveHeaders)
		{
			rec.head = dnaBufferPtr + bufferPos;
			DecompressReadId(rec, mainCtx.id);
			bufferPos += rec.headLen;
		}

		ASSERT(blockDesc.isLenConst);
		rec.seqLen = blockDesc.header.recMinLen;

		ASSERT(bufferPos + rec.seqLen * 2 <= buffer_.size);

		rec.seq = dnaBufferPtr + bufferPos;
		rec.qua = rec.seq + rec.seqLen;


		// decompress dna and qua
		//
		DecompressReadSequence(rec);

		DecompressReadQuality(rec);


		bufferPos += rec.seqLen * 2;			// dna + qua
	}

	ASSERT(bufferPos == buffer_.size);
}


void RawDecompressorSE::DecompressReadSequence(FastqRecord &rec_)
{
	for (uint32 i = 0; i < rec_.seqLen; ++i)
		rec_.seq[i] = mainCtx.dna.rawCoder->GetByte();
}

void RawDecompressorSE::DecompressReadQuality(FastqRecord &rec_)
{
	IQualityStoreBase::DecompressReadQuality(rec_, mainCtx.qua, params);
}


void RawDecompressorSE::DecompressDna(DataChunk& outChunk_, uint64& outSize_, const DataChunk& inChunk_, uint64 inSize_, uint64 inOffset_)
{
	DecompressBuffer(*ppmdCoder, outChunk_, outSize_, inChunk_, inSize_, inOffset_);
}

void RawDecompressorSE::DecompressQuality(DataChunk& outChunk_, uint64& outSize_, const DataChunk& inChunk_, uint64 inSize_, uint64 inOffset_)
{
	// HINT: here can be performed extra quality information compressed
	//
	//. ..
	//


	// compress buffers
	//
	switch (params.quality.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		DecompressBuffer(*ppmdCoder, outChunk_, outSize_, inChunk_, inSize_, inOffset_);
		ASSERT(outSize_ == blockDesc.header.rawDnaStreamSize);
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	case QualityCompressionParams::MET_8BIN:
	{
		// in case of using custom encoder, just copy the decompressed binary data
		// to the output buffer
		std::copy(inChunk_.data.Pointer() + inOffset_, inChunk_.data.Pointer() + inOffset_ + inSize_,
				  outChunk_.data.Pointer());
		outSize_ = inSize_;
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
		// in case of using custom encoder, just copy the decompressed binary data
		// to the output buffer
		std::copy(inChunk_.data.Pointer() + inOffset_, inChunk_.data.Pointer() + inOffset_ + inSize_,
				  outChunk_.data.Pointer());
		outSize_ = inSize_;
		break;
	}
	}
}



LzCompressorPE::LzCompressorPE(const CompressorParams& params_,
							   const QualityCompressionData& globalQuaData_,
							   const FastqRawBlockStats::HeaderStats& headData_,
							   const CompressorAuxParams& auxParams_)
	:	LzCompressorSE(params_, globalQuaData_, headData_, auxParams_)
	,	categorizer(params_.minimizer)
	,	dryFastqBuffer_2(NULL)
	,	dryFastqWriter_2(NULL)
{
	for (uint32 i = 0; i < params_.classifier.maxPairLzWindowSize; ++i)
		pairHistory.push_back(new LzPairMatch());

	if (auxParams_.dry_run)
	{
		dryFastqBuffer_2 = new Buffer(FastqChunk::DefaultBufferSize);
		dryFastqWriter_2 = new BitMemoryWriter(*dryFastqBuffer_2);
	}
}


LzCompressorPE::~LzCompressorPE()
{
	for (auto m : pairHistory)
		delete m;

	if (dryFastqWriter_2 != NULL)
		delete dryFastqWriter_2;

	if (dryFastqBuffer_2 != NULL)
		delete dryFastqBuffer_2;
}

void LzCompressorPE::ClearPairBuffer()
{
	for (auto m : pairHistory)
		m->Clear();

	pairSignatureMap.clear();
}




void LzCompressorPE::SetupBufferMask(std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_)
{
	LzCompressorSE::SetupBufferMask(ppmdBufferCompMask_, buffersNum_);

	ppmdBufferCompMask_[FastqWorkBuffersPE::LetterXBufferPE] = false;
	ppmdBufferCompMask_[FastqWorkBuffersPE::MatchBinaryBufferPE] = false;
	ppmdBufferCompMask_[FastqWorkBuffersPE::FlagBufferPE] = false;
#if (ENC_HR_AC)
	ppmdBufferCompMask_[FastqWorkBuffersPE::HardReadsBufferPE] = false;
#endif
}


void LzCompressorPE::StartEncoding(std::vector<DataChunk*>& buffers_)
{
	LzCompressorSE::StartEncoding(buffers_);

	//pairCtx.lzRle0Coder = new Rle0Encoder(*mainCtx.writers[FastqWorkBuffersPE::LzIdBufferPE]);
	pairCtx.matchRcCoder = new DnaEncoders::RevEncoder(*mainCtx.writers[FastqWorkBuffersPE::MatchBinaryBufferPE]);
	pairCtx.lettersXRcCoder = new DnaEncoders::LettersEncoder(*mainCtx.writers[FastqWorkBuffersPE::LetterXBufferPE]);
	pairCtx.matchRleCoder = new BinaryRleEncoder(*mainCtx.writers[FastqWorkBuffersPE::MatchRLEBufferPE]);
	pairCtx.flagCoder = new DnaEncodersPE::FlagEncoder(*mainCtx.writers[FastqWorkBuffersPE::FlagBufferPE]);
#if (ENC_HR_AC)
	pairCtx.hardCoder = new DnaEncodersPE::HardEncoder(*mainCtx.writers[FastqWorkBuffersPE::HardReadsBufferPE]);
#endif

	// start encoders -- todo: bind rev-rc-coder
	//
	//pairCtx.lzRle0Coder->Start();
	pairCtx.matchRcCoder->Start();
	pairCtx.lettersXRcCoder->Start();
	pairCtx.matchRleCoder->Start();
	pairCtx.flagCoder->Start();
#if (ENC_HR_AC)
	pairCtx.hardCoder->Start();
#endif

	if (auxParams.dry_run)
	{
		dryFastqWriter_2->SetPosition(0);
	}
}


void LzCompressorPE::EndEncoding(std::vector<DataChunk*>& buffers_)
{
	//pairCtx.lzRle0Coder->End();
	pairCtx.matchRcCoder->End();
	pairCtx.lettersXRcCoder->End();
	pairCtx.matchRleCoder->End();
	pairCtx.flagCoder->End();
#if (ENC_HR_AC)
	pairCtx.hardCoder->End();
#endif

	//TFree(pairCtx.lzRle0Coder);
	TFree(pairCtx.matchRcCoder);
	TFree(pairCtx.lettersXRcCoder);
	TFree(pairCtx.matchRleCoder);
	TFree(pairCtx.flagCoder);
#if (ENC_HR_AC)
	TFree(pairCtx.hardCoder);
#endif

	if (auxParams.dry_run)
	{
		// finalize for dry run
		//
		dryFastqWriter_2->Flush();


		// we need mutex here -- to synchronize between paired files
		//
		std::unique_lock<std::mutex> lock(*auxParams.pe_mutex);

		LzCompressorSE::EndEncoding(buffers_);  // <--- writes there dry buffer #1

		// POSIX requires C stdio FILE* ops to be atomic -- no mutex here required
		//
		fwrite(dryFastqWriter_2->Pointer(), 1, dryFastqWriter_2->Position(), auxParams.f_uncompressed_2);
	}
	else
	{
		LzCompressorSE::EndEncoding(buffers_);
	}
}



void LzCompressorPE::CompressMatch(const MatchNode* node_, LzContext& lzEncoder_)
{
	if (auxParams.dry_run)
	{
		const FastqRecord* r = node_->record;
		FastqRecordBuffer buf;
		buf.seqLen = r->seqLen;
		buf.auxLen = r->auxLen;
		buf.headLen = r->headLen;
		buf.head = r->head;

		if (r->IsReadReverse())
		{
			r->ComputeRC(buf);
		}
		else
		{
			buf.CopyFrom(*r);
		}

		if (r->IsPairSwapped())
			buf.SwapReads();

		FastqRecord record = buf;
		FastqRecord pair = record.GetPair();

		std::vector<char> headBuf;
		if (params.archType.readsHaveHeaders && headData.pairedEndFieldIdx != 0)
		{
			// prepare the header
			//
			headBuf.resize(record.headLen);
			pair.head = headBuf.data();

			// find the pair indicator -- copied from FastqParserPE
			//

			// here we set the proper headers
			//
			std::copy(record.head, record.head + record.headLen, pair.head);


			// check whether the /1 /2 fields are present
			//
			if (headData.pairedEndFieldIdx != 0)
			{
				// find position of the /2 token
				const std::string separators = FastqRawBlockStats::HeaderStats::Separators();
				ASSERT(headData.pairedEndFieldIdx != 0);
				uint32 pairTokenPos = 0;

				uint32 fidx = 0;
				for (uint32 i = 0; i <= pair.headLen; ++i)
				{
					if (!std::count(separators.begin(), separators.end(), pair.head[i]) && (i != pair.headLen))
						continue;

					fidx++;

					if (fidx == headData.pairedEndFieldIdx)
					{
						pairTokenPos = i+1;
						break;
					}
				}

				ASSERT(pair.head[pairTokenPos] == '1');
				headBuf[pairTokenPos] = '2';
			}
		}

		// output /1
		//
		record.auxLen = 0;				// we need to reset the aux len to output just 1 record
		DryRunCompress(record, dryFastqWriter);


		// output /2
		//
		DryRunCompress(pair, dryFastqWriter_2);


		// update the stats which are normally handled by LzCompressorSE::CompressRead()
		blockStats->recordsCount++;
	}
	else
	{
		// compress main read -- read ids are processed there
		// -- we use only one and deduce PE from SE read
		//
		LzCompressorSE::CompressMatch(node_, lzEncoder_);


		// compress aux read
		//
		curReadMatchTypeSe = LzMatch;

		CompressPair(*node_->record);
	}


	// update stats
	//
	// ..
}


void LzCompressorPE::CompressRead(const FastqRecord& record_, ReadMatchType type_, bool aux_)
{
	// compress PE
	//
	if (auxParams.dry_run)  // copied from above
	{
		FastqRecordBuffer buf;
		buf.seqLen = record_.seqLen;
		buf.auxLen = record_.auxLen;
		buf.headLen = record_.headLen;
		buf.head = record_.head;

		if (record_.IsReadReverse())
		{
			record_.ComputeRC(buf);
		}
		else
		{
			buf.CopyFrom(record_);
		}

		if (record_.IsPairSwapped())
			buf.SwapReads();

		FastqRecord rec = buf;
		FastqRecord pair = buf.GetPair();

		std::vector<char> headBuf;
		if (params.archType.readsHaveHeaders && headData.pairedEndFieldIdx != 0)
		{
			// prepare the header
			//
			headBuf.resize(rec.headLen);
			pair.head = headBuf.data();

			// find the pair indicator -- copied from FastqParserPE
			//

			// here we set the proper headers
			//
			std::copy(rec.head, rec.head + rec.headLen, pair.head);


			// check whether the /1 /2 fields are present
			//
			if (headData.pairedEndFieldIdx != 0)
			{
				// find position of the /2 token
				const std::string separators = FastqRawBlockStats::HeaderStats::Separators();
				ASSERT(headData.pairedEndFieldIdx != 0);
				uint32 pairTokenPos = 0;

				uint32 fidx = 0;
				for (uint32 i = 0; i <= pair.headLen; ++i)
				{
					if (!std::count(separators.begin(), separators.end(), pair.head[i]) && (i != pair.headLen))
						continue;

					fidx++;

					if (fidx == headData.pairedEndFieldIdx)
					{
						pairTokenPos = i+1;
						break;
					}
				}

				ASSERT(pair.head[pairTokenPos] == '1');
				headBuf[pairTokenPos] = '2';
			}
		}

		// output /1
		//
		rec.auxLen = 0;				// we need to reset the aux len to output just 1 record
		DryRunCompress(rec, dryFastqWriter);


		curReadMatchTypeSe = type_;


		// output /2
		//
		DryRunCompress(pair, dryFastqWriter_2);


		// update the stats which are normally handled by LzCompressorSE::CompressRead()
		blockStats->recordsCount++;
	}
	else
	{
		// compress first read data -- read id will is processed there
		//
		LzCompressorSE::CompressRead(record_, type_, aux_);


		// compress aux read
		//
		curReadMatchTypeSe = type_;

		CompressPair(record_);
	}
}


void LzCompressorPE::CompressPair(const FastqRecord &record_)
{
	// store flags
	//
	// TODO: custom rc
	mainCtx.writers[FastqWorkBuffersPE::SwapFlagBuffer]->PutByte(record_.IsPairSwapped());


	// compute signatures
	//
	FastqRecord pair = record_.GetPair();

	// find signatures from the first and second halves respectively
	//
	auto signatures1 = categorizer.FindMinimizers(pair, 0, pair.seqLen / 2 - (params.minimizer.signatureLen - 1));
	auto signatures2 = categorizer.FindMinimizers(pair, pair.seqLen / 2, 0);
	for (auto kv : signatures1)
	{
		if (signatures2.count(kv.first))
			signatures2.erase(kv.first);
	}

	auto signatures = signatures1;
	signatures.insert(signatures2.begin(), signatures2.end());



	// pop the first element from buffer and remove the signature mapping
	//
	LzPairMatch* lastElem = pairHistory.back();
	pairHistory.pop_back();

	for (uint32 sig : lastElem->signatures)
	{
		if (sig == 0)	// invalid sig
			break;

		ASSERT(pairSignatureMap.count(sig) > 0);

		auto& lzs = pairSignatureMap.at(sig);
		ASSERT(!lzs.empty());

		auto it = std::find(lzs.begin(), lzs.end(), lastElem);
		ASSERT(it != lzs.end());
		lzs.erase(it);			// invalidates elements -- TODO: pointers
	}




	// search for potential lz-matches
	//
	std::tuple<ReadsClassifierSE::MatchResult, const FastqRecord*, uint16> match;
	for (auto& sig : signatures)
	{
		if (pairSignatureMap.count(sig.first) == 0)
			continue;

		// calculate the cost of potential match per each read in the buffer
		//
		const auto& lzs = pairSignatureMap.at(sig.first);
		for (auto lz : lzs)
		{
			FastqRecord plz = lz->rec->GetPair();
			for (uint16 pos : lz->sigPos)
			{
				if (readsClassifier.UpdateLzMatchResult(std::get<0>(match),
														pair.seq, pair.seqLen, sig.second,
														plz.seq, plz.seqLen, pos))
				{
					std::get<1>(match) = lz->rec;
					std::get<2>(match) = pos;
				}
			}
		}
	}

	// lz or hard ?
	const int32 encodeThreshold = (params.classifier.pairEncodeThresholdValue == 0)
				? record_.seqLen / 1.5
				: params.classifier.pairEncodeThresholdValue;

	ReadsClassifierSE::MatchResult& mr = std::get<0>(match);
	const uint32 mism = (mr.cost.cost - ABS(mr.shift * params.classifier.shiftCost)) / params.classifier.mismatchCost;


	uint32 flag = ReadDifficultPE;
	if (mr.cost.cost <= encodeThreshold)		// check better threshold
	{
		// there should not be any exact matches, except duplicates
		ASSERT(mr.cost.cost != 0 || mr.cost.noMismatches);

		// find the read in the buffer and calculate the id
		auto it = std::find_if(pairHistory.begin(), pairHistory.end(),
							   [&match](const LzPairMatch* pm_)
								{return pm_->rec == std::get<1>(match);});
		ASSERT(it != pairHistory.end());

		mr.prevId = it - pairHistory.begin();

		if (mr.cost.noMismatches)
		{
			if (mr.cost.cost == 0)
				flag = ReadIdenticalPE;
			else
				flag = ReadShiftOnlyPE;
		}
		else
		{
			if (mism > params.maxMismatchesLowCost)
				flag = ReadFullExpensivePE;
			else
				flag = ReadFullEncodePE;
		}
	}


	// store flag
	//
	//mainCtx.writers[FastqWorkBuffersPE::FlagBufferPE]->PutByte(flag);
	pairCtx.flagCoder->coder.EncodeSymbol(pairCtx.flagCoder->rc, flag, curReadMatchTypeSe);


	// encode read depending on the flag
	//
	if (flag == ReadDifficultPE)
	{
		// store full sequence
		//
#if (ENC_HR_AC)
		for (int32 i = 0; i < pair.seqLen; ++i)
		{
			uint32 x = dnaToIdx[pair.seq[i]];
			ASSERT(x <= 4);
			pairCtx.hardCoder->coder.EncodeSymbol(pairCtx.hardCoder->rc, x);
		}
#else
		BitMemoryWriter* seqWriter = mainCtx.writers[FastqWorkBuffersPE::HardReadsBufferPE];
		for (int32 i = 0; i < pair.seqLen; ++i)
		{
			seqWriter->PutByte(pair.seq[i]);
		}
#endif
	}
	else if (flag == ReadIdenticalPE)
	{
		// store lz id
		//
		//pairCtx.lzRle0Coder->PutSymbol(std::get<0>(match).prevId);
		mainCtx.writers[FastqWorkBuffersPE::LzIdBufferPE]->Put2Bytes(std::get<0>(match).prevId);
	}
	else
	{
		const ReadsClassifierSE::MatchResult& mr = std::get<0>(match);

		blockStats->freqs["LzMatches_PE-mism"][mism]++;


		// store shift
		//
		ASSERT(ABS(mr.shift) < 127);
		mainCtx.writers[FastqWorkBuffersPE::ShiftBufferPE]->PutByte((int32)(ShiftOffset + mr.shift));

		// store lz id
		//
		//pairCtx.lzRle0Coder->PutSymbol(mr.prevId);
		mainCtx.writers[FastqWorkBuffersPE::LzIdBufferPE]->Put2Bytes(mr.prevId);


		// encode shift differences in respect to lz match
		//
		const LzPairMatch* lz = pairHistory[mr.prevId];
		const FastqRecord lzPair = lz->rec->GetPair();
		const uint16 lzMinPos = std::get<2>(match);

		const char* bestSeq = lzPair.seq;
		uint32 bestLen = lzPair.seqLen;
		uint32 bestPos = lzMinPos;

		const char* newSeq = pair.seq;
		uint32 newLen = pair.seqLen;

		if (mr.shift >= 0)
		{
			ASSERT((int32)bestLen >= mr.shift);
			ASSERT((int32)bestPos >= mr.shift);

			bestSeq += mr.shift;
			bestLen -= mr.shift;
			bestPos -= mr.shift;
		}
		else
		{
			ASSERT(newLen >= newLen);

			// encode insertion
			for (int32 i = 0; i < -mr.shift; ++i)
			{
				pairCtx.lettersXRcCoder->coder.EncodeSymbol(pairCtx.lettersXRcCoder->rc,
															dnaToIdx[(int32)pair.seq[i]], dnaToIdx[(int32)'N']);
			}

			newSeq += -mr.shift;
			newLen -= -mr.shift;

			// the '<=' was changed from '<' due to compression of rev-compl whole trees
			ASSERT((int32)bestPos - mr.shift < (int32)(pair.seqLen));
		}
		uint32 minLen = MIN(bestLen, newLen);


		// encode differences in respect to lz match
		//
		if (flag == ReadFullEncodePE)
		{
			// TODO: optimize this routine
			for (uint32 i = 0; i < minLen; ++i)
			{
				if (bestSeq[i] == newSeq[i])
					pairCtx.matchRleCoder->PutSymbol(true);
				else
				{
					pairCtx.matchRleCoder->PutSymbol(false);
					pairCtx.lettersXRcCoder->coder.EncodeSymbol(pairCtx.lettersXRcCoder->rc,
																dnaToIdx[(int32)newSeq[i]],
																dnaToIdx[(int32)bestSeq[i]]);
				}
			}
		}
		else if (flag == ReadFullExpensivePE)
		{
			for (uint32 i = 0; i < minLen; ++i)
			{
				pairCtx.matchRcCoder->coder.EncodeSymbol(pairCtx.matchRcCoder->rc,
														 bestSeq[i] == newSeq[i]);

				if (bestSeq[i] != newSeq[i])
				{
					pairCtx.lettersXRcCoder->coder.EncodeSymbol(pairCtx.lettersXRcCoder->rc,
																dnaToIdx[(int32)newSeq[i]],
																dnaToIdx[(int32)bestSeq[i]]);
				}
			}
		}


		// encode last shift insertion
		//
		for (uint32 i = minLen; i < newLen; ++i)
		{
			pairCtx.lettersXRcCoder->coder.EncodeSymbol(pairCtx.lettersXRcCoder->rc,
														dnaToIdx[(int32)newSeq[i]],
														dnaToIdx[(int32)'N']);
		}
	}



	// add the read to the history buffer (front) and signature mapping -- TODO: refine this one
	//
	lastElem->Clear();
	lastElem->rec = (FastqRecord*)&record_;


	std::pair<decltype(signatures)*, decltype(signatures)*> ss = (signatures1.size() > signatures2.size()) ?
				std::make_pair(&signatures2, &signatures1) :
				std::make_pair(&signatures1, &signatures2);

	std::vector<decltype(signatures)::value_type> fs;
	uint32 numSig = 0;
	decltype(signatures)::const_iterator si = ss.first->begin();
	for ( ; numSig < MIN(2, ss.first->size()); numSig++, si++ )
		fs.push_back(*si);

	si = ss.second->begin();
	for ( ; numSig < MIN(4, ss.first->size() + ss.second->size()); numSig++, si++ )
		fs.push_back(*si);

	for (uint32 i = 0; i < fs.size(); ++i)
	{
		ASSERT(fs[i].first != 0);

		lastElem->signatures[i] = fs[i].first;
		lastElem->sigPos[i] = fs[i].second;

		// add pointer to the new element for reverse search
		pairSignatureMap[fs[i].first].push_back(lastElem);
	}


	// update the prev-buffer
	//
	if (flag != ReadIdenticalPE)
		pairHistory.push_front(lastElem);
	else
		pairHistory.push_back(lastElem);


	// compress quality
	//
	pair.SetReadReverse(record_.IsReadReverse());
	CompressReadQuality(pair);



	// update stats
	//
	switch (flag)
	{
	case ReadDifficultPE:
		blockStats->counts["HardReads_PE"]++;
		break;

	case ReadIdenticalPE:
		blockStats->counts["ExactMatches_PE"]++;
		break;

	case ReadFullEncodePE:
	case ReadFullExpensivePE:
	case ReadShiftOnlyPE:
	{
		const ReadsClassifierSE::MatchResult& mr = std::get<0>(match);

		blockStats->counts["LzMatches_PE"]++;
		blockStats->freqs["LzId_PE"][mr.prevId]++;
		blockStats->freqs["Shift_PE"][mr.shift]++;
		blockStats->freqs["Cost_PE"][mr.cost.cost]++;
		break;
	}
	}

	switch (curReadMatchTypeSe)
	{
	case HardRead:
		blockStats->freqs["SE_HR-2-PE-flag"][flag]++;
		break;

	case ExactMatch:
		blockStats->freqs["SE_EM-2-PE-flag"][flag]++;
		break;

	case LzMatch:
		blockStats->freqs["SE_LZ-2-PE-flag"][flag]++;
		break;

	case ContigRead:
		blockStats->freqs["SE_CR-2-PE-flag"][flag]++;
		break;
	}
}



void LzDecompressorPE::Decompress(CompressedFastqBlock& compBin_,
								  std::vector<FastqRecord>& reads_,
								  FastqCompressedBin& fastqWorkBin_,
								  FastqChunk& dnaBuffer_)
{
	pairHistory.clear();

#if EXTRA_MEM_OPT
	pairHistory.clear();
#endif

	LzDecompressorSE::Decompress(compBin_, reads_, fastqWorkBin_, dnaBuffer_);
}


void LzDecompressorPE::DecompressMetaAndLinkRecords(std::vector<FastqRecord> &reads_, FastqChunk &chunk_)
{
	// TODO: now only supports const-length records
	//
	ASSERT(blockDesc.header.recMinLen == blockDesc.header.recMaxLen);

	char* dnaBufferPtr = (char*)chunk_.data.Pointer();
	uint64 pos = 0;

	for (uint64 i = 0; i < reads_.size(); ++i)
	{
		FastqRecord& r = reads_[i];

		// decompress headers
		//
		if (params.archType.readsHaveHeaders)
		{
			r.head = dnaBufferPtr + pos;
			DecompressReadId(r, mainCtx.id);
			pos += r.headLen;
		}

		// decompress read length
		//
		r.seqLen = r.auxLen = blockDesc.header.recMaxLen;
		r.seq = dnaBufferPtr + pos;
		r.qua = r.seq + r.seqLen + r.auxLen;

		pos += (r.seqLen + r.auxLen) * 2;
	}

	chunk_.size = pos;
}


void LzDecompressorPE::SetupBufferMask(std::vector<bool>& ppmdBufferCompMask_, uint32 buffersNum_)
{
	LzDecompressorSE::SetupBufferMask(ppmdBufferCompMask_, buffersNum_);

	ppmdBufferCompMask_[FastqWorkBuffersPE::LetterXBufferPE] = false;
	ppmdBufferCompMask_[FastqWorkBuffersPE::MatchBinaryBufferPE] = false;
	ppmdBufferCompMask_[FastqWorkBuffersPE::FlagBufferPE] = false;
#if (ENC_HR_AC)
	ppmdBufferCompMask_[FastqWorkBuffersPE::HardReadsBufferPE] = false;
#endif
}


void LzDecompressorPE::StartDecoding(std::vector<DataChunk*>& buffers_)
{
	LzDecompressorSE::StartDecoding(buffers_);

	// create decoders -- TODO: create only once and reuse
	//
	//pairCtx.lzRle0Coder = new Rle0Decoder(*mainCtx.readers[FastqWorkBuffersPE::LzIdBufferPE]);
	pairCtx.matchRcCoder = new DnaDecoders::RevDecoder(*mainCtx.readers[FastqWorkBuffersPE::MatchBinaryBufferPE]);
	pairCtx.lettersXRcCoder = new DnaDecoders::LettersDecoder(*mainCtx.readers[FastqWorkBuffersPE::LetterXBufferPE]);
	pairCtx.matchRleCoder = new BinaryRleDecoder(*mainCtx.readers[FastqWorkBuffersPE::MatchRLEBufferPE]);
	pairCtx.flagCoder = new DnaDecodersPE::FlagDecoder(*mainCtx.readers[FastqWorkBuffersPE::FlagBufferPE]);
#if (ENC_HR_AC)
	pairCtx.hardCoder = new DnaDecodersPE::HardDecoder(*mainCtx.readers[FastqWorkBuffersPE::HardReadsBufferPE]);
#endif

	// start decoders
	//
	//pairCtx.lzRle0Coder->Start();
	pairCtx.matchRcCoder->Start();
	pairCtx.lettersXRcCoder->Start();
	pairCtx.matchRleCoder->Start();
	pairCtx.flagCoder->Start();
#if (ENC_HR_AC)
	pairCtx.hardCoder->Start();
#endif
}


void LzDecompressorPE::EndDecoding()
{
	//pairCtx.lzRle0Coder->End();
	pairCtx.matchRcCoder->End();
	pairCtx.lettersXRcCoder->End();
	pairCtx.matchRleCoder->End();
	pairCtx.flagCoder->End();
#if (ENC_HR_AC)
	pairCtx.hardCoder->End();
#endif

	//TFree(pairCtx.lzRle0Coder);
	TFree(pairCtx.matchRcCoder);
	TFree(pairCtx.lettersXRcCoder);
	TFree(pairCtx.matchRleCoder);
	TFree(pairCtx.flagCoder);
#if (ENC_HR_AC)
	TFree(pairCtx.hardCoder);
#endif

	LzDecompressorSE::EndDecoding();
}


void LzDecompressorPE::DecompressRecord(FastqRecord& rec_, ReadFlags flag_,
										const FastqRecord* auxRec_, bool aux_)
{
	LzDecompressorSE::DecompressRecord(rec_, flag_, auxRec_, aux_);

	// TODO: LUT value mapping
	//
	switch (flag_)
	{
	case ReadDifficult:
		curReadMatchTypeSe = HardRead;
		break;

	case ReadIdentical:
		curReadMatchTypeSe = ExactMatch;
		break;

	case ReadShiftOnly:
	case ReadFullEncode:
	case ReadFullExpensive:
		curReadMatchTypeSe = LzMatch;
		break;

	case ReadContigGroupNext:
		curReadMatchTypeSe = ContigRead;
		break;

	default:
		ASSERT(0);
		break;
	}


	// decompress PE
	//
	DecompressPair(rec_);
}


void LzDecompressorPE::DecompressPair(FastqRecord &record_)
{
	// read swap flag
	//
	record_.SetPairSwapped(mainCtx.readers[FastqWorkBuffersPE::SwapFlagBuffer]->GetByte() != 0);


	// read encoding flag
	//
	//uint32 flag = mainCtx.readers[FastqWorkBuffersPE::FlagBufferPE]->GetByte();
	uint32 flag = pairCtx.flagCoder->coder.DecodeSymbol(pairCtx.flagCoder->rc, curReadMatchTypeSe);

	ASSERT(flag == ReadShiftOnlyPE || flag == ReadFullEncodePE || flag == ReadFullExpensivePE
		   || flag == ReadDifficultPE || flag == ReadIdenticalPE);



	FastqRecord pair = record_.GetPair();

	if (flag == ReadDifficultPE)
	{
		// read whole sequence
		//
#if (ENC_HR_AC)
		for (int32 i = 0; i < pair.seqLen; ++i)
		{
			int32 x = pairCtx.hardCoder->coder.DecodeSymbol(pairCtx.hardCoder->rc);
			ASSERT(x <= 4);
			char c = idxToDna[x];
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			pair.seq[i] = (char)c;
		}
#else
		BitMemoryReader* seqReader = mainCtx.readers[FastqWorkBuffersPE::HardReadsBufferPE];

		for (int32 i = 0; i < pair.seqLen; ++i)
		{
			int32 c = seqReader->GetByte();
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			pair.seq[i] = (char)c;
		}
#endif
	}
	else if (flag == ReadIdenticalPE)
	{
		//int32 prevId = pairCtx.lzRle0Coder->GetSym();
		int32 prevId = mainCtx.readers[FastqWorkBuffersPE::LzIdBufferPE]->Get2Bytes();
		ASSERT(prevId < (int32)pairHistory.size());

		const FastqRecord* lzRec = (const FastqRecord*)*(pairHistory.rbegin() + prevId);
		const FastqRecord lzPair = lzRec->GetPair();

		std::copy(lzPair.seq, lzPair.seq + lzPair.seqLen, pair.seq);
	}
	else
	{
		int32 recLen = pair.seqLen;
		char* recStr = pair.seq;

		int32 shift = (int32)(mainCtx.readers[FastqWorkBuffersPE::ShiftBufferPE]->GetByte()) - ShiftOffset;

		ASSERT(ABS(shift) < recLen - params.minimizer.signatureLen);

		//int32 prevId = pairCtx.lzRle0Coder->GetSym();
		int32 prevId = mainCtx.readers[FastqWorkBuffersPE::LzIdBufferPE]->Get2Bytes();
		ASSERT(prevId < (int32)pairHistory.size());

		const FastqRecord* lzRec = (const FastqRecord*)*(pairHistory.rbegin() + prevId);
		const FastqRecord lzPair = lzRec->GetPair();

		const char* bestSeq = lzPair.seq;
		int32 bestLen = lzPair.seqLen;
		//int32 bestPos = lzPair.minimPos;

		if (shift >= 0)
		{
			bestSeq += shift;
			bestLen -= shift;
			//bestPos -= shift;
		}
		else
		{
			for (int32 i = 0; i < -shift; ++i)
			{
				int32 c = pairCtx.lettersXRcCoder->coder.DecodeSymbol(pairCtx.lettersXRcCoder->rc,
																	  dnaToIdx[(int)'N']);
				ASSERT(c < 5);
				c = idxToDna[c];

				ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
				pair.seq[i] = c;
			}

			recStr += -shift;
			recLen -= -shift;
		}

		int32 minLen = MIN(bestLen, recLen);
		ASSERT(minLen > 0);


		if (flag == ReadFullEncodePE)
		{
			// TODO: optimize this routine
			for (int32 i = 0; i < minLen; ++i)
			{
				if (pairCtx.matchRleCoder->GetSym())
				{
					recStr[i] = bestSeq[i];
				}
				else
				{
					int32 c = pairCtx.lettersXRcCoder->coder.DecodeSymbol(pairCtx.lettersXRcCoder->rc,
																		  dnaToIdx[(int32)bestSeq[i]]);
					ASSERT(c < 5);
					c = idxToDna[c];

					ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
					recStr[i] = c;
				}
			}
		}
		else if (flag == ReadFullExpensivePE)
		{
			for (int32 i = 0; i < minLen; ++i)
			{
				if (pairCtx.matchRcCoder->coder.DecodeSymbol(pairCtx.matchRcCoder->rc) != 0)
				{
					recStr[i] = bestSeq[i];
				}
				else
				{
					int32 c = pairCtx.lettersXRcCoder->coder.DecodeSymbol(pairCtx.lettersXRcCoder->rc,
																		  dnaToIdx[(int32)bestSeq[i]]);
					ASSERT(c < 5);
					c = idxToDna[c];

					ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
					recStr[i] = c;
				}
			}
		}
		else					// all symbols match
		{
			// TODO: use std::copy
			for (int32 i = 0; i < minLen; ++i)
				recStr[i] = bestSeq[i];
		}


		// read the remaining part
		//
		for (int32 i = minLen; i < recLen; ++i)
		{
			int32 c = pairCtx.lettersXRcCoder->coder.DecodeSymbol(pairCtx.lettersXRcCoder->rc,
																  dnaToIdx[(int)'N']);
			ASSERT(c < 5);
			c = idxToDna[c];

			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			recStr[i] = c;
		}
	}


	// update LZ context
	if (flag != ReadIdenticalPE)
		pairHistory.push_back((const FastqRecord*)&record_);


	// decompress quality
	//
	pair.SetReadReverse(record_.IsReadReverse());
	DecompressReadQuality(pair);
}



RawCompressorPE::RawCompressorPE(const CompressorParams &params_,
								 const QualityCompressionData &globalQuaData_,
								 const FastqRawBlockStats::HeaderStats &headData_,
								 const CompressorAuxParams &auxParams_)
	:	RawCompressorSE(params_, globalQuaData_, headData_, auxParams_)
	,	dryFastqBuffer_2(NULL)
	,	dryFastqWriter_2(NULL)
{
	if (auxParams_.dry_run)
	{
		dryFastqBuffer_2 = new Buffer(FastqChunk::DefaultBufferSize);
		dryFastqWriter_2 = new BitMemoryWriter(*dryFastqBuffer_2);
	}
}

RawCompressorPE::~RawCompressorPE()
{
	if (dryFastqWriter_2 != NULL)
		delete dryFastqWriter_2;

	if (dryFastqBuffer_2 != NULL)
		delete dryFastqBuffer_2;
}


void RawCompressorPE::DryRunCompress(const FastqRecord& record_, BitMemoryWriter* )
{
	// copied from LzCompressor -- TODO: merge into one method
	//

	FastqRecord rec = record_;
	FastqRecord pair = rec.GetPair();

	pair.SetReadReverse(rec.IsReadReverse());

	if (rec.IsPairSwapped())
	{
		std::swap(rec, pair);
	}

	std::vector<char> headBuf;
	if (params.archType.readsHaveHeaders && headData.pairedEndFieldIdx != 0)
	{
		// prepare the header
		//
		headBuf.resize(rec.headLen);
		pair.head = headBuf.data();

		// find the pair indicator -- copied from FastqParserPE
		//

		// here we set the proper headers
		//
		std::copy(rec.head, rec.head + rec.headLen, pair.head);


		// check whether the /1 /2 fields are present
		//
		if (headData.pairedEndFieldIdx != 0)
		{
			// find position of the /2 token
			const std::string separators = FastqRawBlockStats::HeaderStats::Separators();
			ASSERT(headData.pairedEndFieldIdx != 0);
			uint32 pairTokenPos = 0;

			uint32 fidx = 0;
			for (uint32 i = 0; i <= pair.headLen; ++i)
			{
				if (!std::count(separators.begin(), separators.end(), pair.head[i]) && (i != pair.headLen))
					continue;

				fidx++;

				if (fidx == headData.pairedEndFieldIdx)
				{
					pairTokenPos = i+1;
					break;
				}
			}

			ASSERT(pair.head[pairTokenPos] == '1');
			headBuf[pairTokenPos] = '2';
		}
	}

	// output /1
	//
	rec.auxLen = 0;
	RawCompressorSE::DryRunCompress(rec, dryFastqWriter);


	// output /2
	//
	RawCompressorSE::DryRunCompress(pair, dryFastqWriter_2);
}

void RawCompressorPE::StartEncoding(std::vector<DataChunk *> &buffers_)
{
	RawCompressorSE::StartEncoding(buffers_);

	if (auxParams.dry_run)
	{
		dryFastqWriter_2->SetPosition(0);
	}
}

void RawCompressorPE::EndEncoding(std::vector<DataChunk *> &buffers_)
{
	if (auxParams.dry_run)
	{
		// finalize for dry run
		//
		dryFastqWriter_2->Flush();


		// we need mutex here -- to synchronize between paired files
		//
		std::unique_lock<std::mutex> lock(*auxParams.pe_mutex);
		RawCompressorSE::EndEncoding(buffers_);  // <--- writes there dry buffer #1

		// POSIX requires C stdio FILE* ops to be atomic -- no mutex here required
		//
		fwrite(dryFastqWriter_2->Pointer(), 1, dryFastqWriter_2->Position(), auxParams.f_uncompressed_2);
	}
	else
	{
		RawCompressorSE::EndEncoding(buffers_);
	}
}


void RawCompressorPE::CompressReadSequence(const FastqRecord &rec_)
{
	FastqRecord pair = rec_.GetPair();
	RawCompressorSE::CompressReadSequence(rec_);
	RawCompressorSE::CompressReadSequence(pair);
}

void RawCompressorPE::CompressReadQuality(const FastqRecord &rec_)
{
	RawCompressorSE::CompressReadQuality(rec_);

	FastqRecord pair = rec_.GetPair();
	// read direction is used when encoding/decoding qualities
	pair.SetReadReverse(rec_.IsReadReverse());
	RawCompressorSE::CompressReadQuality(pair);
}



void RawDecompressorPE::DecompressRecords(std::vector<FastqRecord>& reads_, FastqChunk& buffer_)
{
	ASSERT(blockDesc.isLenConst);

	char* dnaBufferPtr = (char*)buffer_.data.Pointer();
	uint64 bufferPos = 0;

	for (FastqRecord& rec : reads_)
	{
		rec.Reset();

		// decompress header
		//
		if (params.archType.readsHaveHeaders)
		{
			rec.head = dnaBufferPtr + bufferPos;
			DecompressReadId(rec, mainCtx.id);
			bufferPos += rec.headLen;
		}

		if (blockDesc.isLenConst)
		{
			rec.seqLen = rec.auxLen = blockDesc.header.recMinLen;
		}
		ASSERT(bufferPos + (rec.seqLen + rec.auxLen * 2) <= buffer_.size);


		// link dna and quality to output chunk buffer
		//
		rec.seq = dnaBufferPtr + bufferPos;
		rec.qua = rec.seq + rec.seqLen + rec.auxLen ;


		// decompress reads
		//
		DecompressReadSequence(rec);

		DecompressReadQuality(rec);

		bufferPos += (rec.seqLen + rec.auxLen) * 2;
	}

	ASSERT(bufferPos == buffer_.size);
}


void RawDecompressorPE::DecompressReadSequence(FastqRecord &rec_)
{
	FastqRecord pair = rec_.GetPair();

	RawDecompressorSE::DecompressReadSequence(rec_);
	RawDecompressorSE::DecompressReadSequence(pair);
}

void RawDecompressorPE::DecompressReadQuality(FastqRecord &rec_)
{
	RawDecompressorSE::DecompressReadQuality(rec_);

	FastqRecord pair = rec_.GetPair();
	// read direction is used when encoding/decoding qualities
	pair.SetReadReverse(rec_.IsReadReverse());
	RawDecompressorSE::DecompressReadQuality(pair);
}



// dynamic decompressor
//
#if 0
void LzDecompressorDynSE::StartDecompress(CompressedFastqBlock &compBin_, FastqCompressedBin &fastqWorkBin_)
{

	// read header and footer
	//
	const uint32 buffersNum = fastqWorkBin_.buffers.size();
	blockDesc.Reset(buffersNum);
	ReadHeader(blockDesc.header, compBin_.dataBuffer);
	blockDesc.isLenConst = (blockDesc.header.recMinLen == blockDesc.header.recMaxLen);

	{
		BitMemoryReader reader(compBin_.dataBuffer.data,
							   blockDesc.header.footerOffset + blockDesc.header.footerSize,
							   blockDesc.header.footerOffset);
		ReadRawFooter(blockDesc.footer, reader);
	}


	// setup and decompress buffers
	//
	SetupBuffers(dnaBuffer_.data, mainCtx.bufferCompMask, buffersNum);


	DecompressBuffers(fastqWorkBin_.buffers, mainCtx.bufferCompMask,
					  compBin_.dataBuffer, LzBlockHeader::Size(buffersNum));

	// start decoding
	//
	StartDecoding(fastqWorkBin_.buffers);

	decompressedBytes = 0;
	decompressedRecs = 0;

	/*

	// -- decompress here the headers -- like in DSRC ;---)
	DecompressMetaAndLinkRecords(reads_, dnaBuffer_);

	//dnaBuffer_.size = blockDesc.header.rawDnaStreamSize * 2;


	// decode records
	//
	DecompressRecords(reads_);


	// finish decoding
	//
	EndDecoding();

	ASSERT(dnaBuffer_.size == blockDesc.header.rawDnaStreamSize * 2 + blockDesc.header.rawIdStreamSize);
	*/
}

void LzDecompressorDynSE::SetupBuffers(Buffer &outBuffer_, std::vector<bool> &ppmdBufferCompMask_, uint32 buffersNum_)
{
	SetupBufferMask(ppmdBufferCompMask_, buffersNum_);
}

void LzDecompressorDynSE::FinishDecompress()
{
	EndDecoding();

	ASSERT(decompressedBytes == blockDesc.header.rawDnaStreamSize * 2 + blockDesc.header.rawIdStreamSize);
}


uint64 LzDecompressorDynSE::DecompressNext(std::vector<FastqRecord> &reads_, FastqChunk &chunk_)
{

	// TODO: now only supports const-length records
	//
	ASSERT(blockDesc.header.recMinLen == blockDesc.header.recMaxLen);

	// TODO: decompress read lengths here
	//

	reads_.clear();

	char* dnaBufferPtr = (char*)chunk_.data.Pointer();
	uint64 pos = 0;

	while (pos + 1024 < chunk_.data.Size())
	{
		FastqRecord r;
		r.seqLen = blockDesc.header.recMaxLen;


		// decompress header and link record with buffers
		//
		if (params.archType.readsHaveHeaders)
		{
			r.head = dnaBufferPtr + pos;
			DecompressReadId(r, mainCtx.id);
			pos += r.headLen;
		}

		// decompress read length
		//
		r.seq = dnaBufferPtr + pos;
		r.qua = r.seq + r.seqLen;

		pos += r.seqLen * 2;


		// decompress record
		//
		ASSERT(consDecoders.size() == 0);

		uint32 flag = mainCtx.readers[FastqWorkBuffersSE::FlagBuffer]->GetByte();
		ASSERT(flag <= ReadContigGroupStart);

		// INFO: we are incrementing the read idx only when decompressing normal reads
		// the groups: sub-tree and contig groups handle the index incrementations
		// by themselves
		// TODO: create the records iterator and pass it
		if (flag == ReadContigGroupStart)
		{
			DecompressConsensus(reads_, recIdx);
		}
		else if (flag == ReadTreeGroupStart)
		{
			DecompressTree(reads_[recIdx - 1],
						   reads_,
						   recIdx);
		}
		else if (flag == ReadIdentical)
		{
			ASSERT(lzContextStack.size() > 0);
			LzContext& lzCtx = lzContextStack.top();

			const FastqRecord* match = *(lzCtx.history.rbegin());

			DecompressRecord(reads_[recIdx], (ReadFlags)flag, match);
			recIdx++;
		}
		else
		{
			ASSERT(flag < ReadTreeGroupStart);
			DecompressRecord(reads_[recIdx], (ReadFlags)flag);
			recIdx++;
		}

	}

	chunk_.size = pos;

}


#endif



// DNA compressor proxy
//

FastqCompressor::FastqCompressor(const CompressorParams& params_,
								 const QualityCompressionData& globalQuaData_,
								 const FastqRawBlockStats::HeaderStats& headData_,
								 const CompressorAuxParams& auxParams_)
	:	params(params_)
	,	globalQuaData(globalQuaData_)
	,	headData(headData_)
	,	auxParams(auxParams_)
	,	lzCompressor(NULL)
	,	rawCompressor(NULL)
{}

FastqCompressor::~FastqCompressor()
{
	TFree(lzCompressor);
	TFree(rawCompressor);
}

void FastqCompressor::Compress(const std::vector<FastqRecord>& reads_,
								 PackContext& packCtx_,
								 uint32 minimizerId_,
								 uint64 rawDnaStreamSize_,
								 FastqCompressedBin& fastqWorkBin_,
								 CompressedFastqBlock &compBin_)
{
	if (minimizerId_ != params.minimizer.SignatureN())
	{
		if (lzCompressor == NULL)
		{
			lzCompressor = params.archType.readType == ArchiveType::READ_SE
					? new LzCompressorSE(params, globalQuaData, headData, auxParams)
					: new LzCompressorPE(params, globalQuaData, headData, auxParams);
		}
		lzCompressor->Compress(reads_, packCtx_, minimizerId_,
							   rawDnaStreamSize_, fastqWorkBin_, compBin_);
	}
	else
	{
		if (rawCompressor == NULL)
		{
			rawCompressor = params.archType.readType == ArchiveType::READ_SE
					? new RawCompressorSE(params, globalQuaData, headData, auxParams)
					: new RawCompressorPE(params, globalQuaData, headData, auxParams);
		}
		rawCompressor->Compress(reads_, packCtx_, minimizerId_,
								rawDnaStreamSize_, fastqWorkBin_, compBin_);
	}
}





FastqDecompressor::FastqDecompressor(const CompressorParams& params_,
								 const QualityCompressionData& globalQuaData_,
								 const FastqRawBlockStats::HeaderStats& headData_,
								 const CompressorAuxParams& auxParams_)
	:	params(params_)
	,	globalQuaData(globalQuaData_)
	,	headData(headData_)
	,	auxParams(auxParams_)
	,	lzDecompressor(NULL)
	,	rawDecompressor(NULL)
{}

FastqDecompressor::~FastqDecompressor()
{
	TFree(lzDecompressor);
	TFree(rawDecompressor);
}



void FastqDecompressor::Decompress(CompressedFastqBlock& compBin_,
								   std::vector<FastqRecord>& reads_,
								   FastqCompressedBin& fastqWorkBin_,
								   FastqChunk& dnaBuffer_)
{
	if (compBin_.signatureId != params.minimizer.SignatureN())
	{
		if (lzDecompressor == NULL)
		{
			lzDecompressor = params.archType.readType == ArchiveType::READ_SE
					? new LzDecompressorSE(params, globalQuaData, headData)
					: new LzDecompressorPE(params, globalQuaData, headData);
		}
		lzDecompressor->Decompress(compBin_, reads_, fastqWorkBin_, dnaBuffer_);
	}
	else
	{
		if (rawDecompressor == NULL)
		{
			rawDecompressor = params.archType.readType == ArchiveType::READ_SE
					? new RawDecompressorSE(params, globalQuaData, headData)
					: new RawDecompressorPE(params, globalQuaData, headData);
		}
		rawDecompressor->Decompress(compBin_, reads_, fastqWorkBin_, dnaBuffer_);
	}
}

