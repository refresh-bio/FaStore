/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"
#include "FastqPacker.h"
#include "BitMemory.h"
#include "BinBlockData.h"
#include "Utils.h"

#include <algorithm>


IFastqPacker::IFastqPacker(const BinModuleConfig& binConfig_)
	:	binConfig(binConfig_)
{
	// dna translation tables
	//
	std::fill(dnaToIdx.begin(), dnaToIdx.end(), -1);
	for (uint32 i = 0; i < 5; ++i)
	{
		int32 c = binConfig.minimizer.dnaSymbolOrder[i];
		ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
		dnaToIdx[c] = i;
	}

	std::fill(idxToDna.begin(), idxToDna.end(), -1);
	for (uint32 i = 0; i < 5; ++i)
	{
		idxToDna[i] = binConfig.minimizer.dnaSymbolOrder[i];
	}


	// qualities translation tables -- those should be used in separate modules
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
}


bool IFastqPacker::ReadNextRecord(BitMemoryReader& metaReader_,
									 BitMemoryReader& dnaReader_,
									 BitMemoryReader& quaReader_,
									 BitMemoryReader& headReader_,
									 const BinPackSettings& settings_,
									 FastqRecord& rec_)
{
	if (dnaReader_.Position() >= dnaReader_.Size())
		return false;

	// read general record info
	//
	const bool hasMinimizer = settings_.suffixLen != 0;
	if (hasMinimizer)
	{
		rec_.SetReadReverse(metaReader_.GetBit() != 0);
		rec_.minimPos = metaReader_.GetBits(LenBits);
	}
	else
	{
		rec_.SetReadReverse(false);			// TODO: check and remove
		rec_.minimPos = 0;
	}

	// read sequence
	//
	ReadDna(metaReader_, dnaReader_, settings_, rec_);


	// read quality
	//
	ReadQuality(metaReader_, quaReader_, settings_, rec_);


	// read header
	//
	if (settings_.usesHeaders)
	{
		ReadHeader(metaReader_, headReader_, settings_, rec_);
	}

	return true;
}


void IFastqPacker::StoreNextRecord(BitMemoryWriter& metaWriter_,
									  BitMemoryWriter& dnaWriter_,
									  BitMemoryWriter& quaWriter_,
									  BitMemoryWriter& headWriter_,
									  const BinPackSettings& settings_,
									  const FastqRecord& rec_)
{
	ASSERT(rec_.seqLen > 0);

	// store meta info
	//
	const bool hasMinimizer = settings_.suffixLen != 0;
	if (hasMinimizer)
	{
		metaWriter_.PutBit(rec_.IsReadReverse());
		metaWriter_.PutBits(rec_.minimPos, LenBits);
	}
	else
	{
		ASSERT(!rec_.IsReadReverse());
		ASSERT(rec_.minimPos == 0);
	}


	// store sequence
	//
	StoreDna(metaWriter_, dnaWriter_, settings_, rec_);


	// store quality
	//
	StoreQuality(metaWriter_, quaWriter_, settings_, rec_);


	// store header
	//
	if (settings_.usesHeaders)
	{
		StoreHeader(metaWriter_, headWriter_, settings_, rec_);
	}
}



void IFastqPacker::StoreDna(BitMemoryWriter& metaWriter_,
								  BitMemoryWriter& dnaWriter_,
							   const BinPackSettings& settings_,
								  const FastqRecord& rec_)
{
	// store sequence
	//
	const bool isDnaPlain = std::find(rec_.seq, rec_.seq + rec_.seqLen, 'N') == rec_.seq + rec_.seqLen;
	metaWriter_.PutBit(isDnaPlain);


	// when saving record, skip writing bytes identifying the signature
	//
	if (isDnaPlain)
	{
		for (uint32 i = 0 ; i < rec_.minimPos; ++i)
		{
			char c = rec_.seq[i];
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T');
			dnaWriter_.Put2Bits(dnaToIdx[(uint32)c]);
		}

		for (uint32 i = rec_.minimPos + settings_.suffixLen; i < rec_.seqLen; ++i)
		{
			char c = rec_.seq[i];
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T');
			dnaWriter_.Put2Bits(dnaToIdx[(uint32)c]);
		}
	}
	else
	{
		for (uint32 i = 0; i < rec_.minimPos; ++i)
		{
			char c = rec_.seq[i];
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			dnaWriter_.PutBits(dnaToIdx[(uint32)c], 3);
		}

		for (uint32 i = rec_.minimPos + settings_.suffixLen; i < rec_.seqLen; ++i)
		{
			char c = rec_.seq[i];
			ASSERT(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
			dnaWriter_.PutBits(dnaToIdx[(uint32)c], 3);
		}
	}
}


void IFastqPacker::StoreQuality(BitMemoryWriter& /*metaWriter_*/,
									  BitMemoryWriter& quaWriter_,
								   const BinPackSettings& /*settings_*/,
									  const FastqRecord& rec_)
{
	// TODO: create a separate module with quality packing
	//
	ASSERT(rec_.qua != NULL);

	const uint32 qbits = binConfig.quaParams.BitsPerBase();
	switch (binConfig.quaParams.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 c = rec_.qua[i];
			ASSERT(c >= binConfig.archiveType.qualityOffset && c < binConfig.archiveType.qualityOffset + 64U);
			c -= binConfig.archiveType.qualityOffset;
			quaWriter_.PutBits(c, qbits);
		}
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		ASSERT(binConfig.quaParams.binaryThreshold < 64);
		ASSERT(qbits == 1);

		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 c = rec_.qua[i];
			ASSERT(c >= binConfig.archiveType.qualityOffset && c < binConfig.archiveType.qualityOffset + 64U);
			c -= binConfig.archiveType.qualityOffset;
			quaWriter_.PutBit(c >= binConfig.quaParams.binaryThreshold);
		}
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 c = rec_.qua[i];
			ASSERT(c >= binConfig.archiveType.qualityOffset && c < binConfig.archiveType.qualityOffset + 64U);
			c = quaToIdx_8bin[c-binConfig.archiveType.qualityOffset];
			quaWriter_.PutBits(c, qbits);
		}
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
        for (uint32 i = 0; i < rec_.seqLen; ++i)
        {
            uint32 c = rec_.qua[i];
			ASSERT(c >= binConfig.archiveType.qualityOffset && c < binConfig.archiveType.qualityOffset + 64U);
			c -= binConfig.archiveType.qualityOffset;
            quaWriter_.PutBits(c, qbits);
        }
        break;
	}

	}
}


void IFastqPacker::StoreHeader(BitMemoryWriter& /*metaWriter_*/,
							   BitMemoryWriter& headWriter_,
							   const BinPackSettings& /*settings_*/,
							   const FastqRecord& rec_)
{
	// now just store the header as it is -- TODO: more advanced compression
	ASSERT(rec_.head != NULL);
	ASSERT(rec_.headLen > 0);
	ASSERT(rec_.head[0] == '@');

	headWriter_.PutBits(rec_.headLen, 8);			// max 255 length
	for (uint32 i = 1; i < rec_.headLen; ++i)
	{
		headWriter_.PutBits(rec_.head[i], 7);		// max 128 ASCII values
	}
}


void IFastqPacker::ReadDna(BitMemoryReader& metaReader_,
								 BitMemoryReader& dnaReader_,
							  const BinPackSettings& settings_,
								 FastqRecord& rec_)
{
	const bool isDnaPlain = metaReader_.GetBit() != 0;

	if (isDnaPlain)
	{
		for (uint32 i = 0; i < rec_.minimPos; ++i)
		{
			rec_.seq[i] = idxToDna[dnaReader_.Get2Bits()];
			ASSERT(rec_.seq[i] != -1);
		}

		for (uint32 i = rec_.minimPos + settings_.suffixLen; i < rec_.seqLen; ++i)
		{
			rec_.seq[i] = idxToDna[dnaReader_.Get2Bits()];
			ASSERT(rec_.seq[i] != -1);
		}
	}
	else
	{
		for (uint32 i = 0; i < rec_.minimPos; ++i)
		{
			uint32 c = dnaReader_.GetBits(3);
			ASSERT(c < 5);
			rec_.seq[i] = idxToDna[c];
			ASSERT(rec_.seq[i] != -1);
		}

		for (uint32 i = rec_.minimPos + settings_.suffixLen; i < rec_.seqLen; ++i)
		{
			uint32 c = dnaReader_.GetBits(3);
			ASSERT(c < 5);
			rec_.seq[i] = idxToDna[c];
			ASSERT(rec_.seq[i] != -1);
		}
	}
}

void IFastqPacker::ReadQuality(BitMemoryReader& /*metaReader_*/,
									 BitMemoryReader& quaReader_,
								  const BinPackSettings& /*settings_*/,
									 FastqRecord& rec_)
{

	// select the apprpriate scheme to read quality
	//
	ASSERT(rec_.qua != NULL);

	const uint32 qbits = binConfig.quaParams.BitsPerBase();
	switch(binConfig.quaParams.method)
	{
	case QualityCompressionParams::MET_NONE:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 c = quaReader_.GetBits(qbits);
			rec_.qua[i] = (char)(c + binConfig.archiveType.qualityOffset);
		}
		break;
	}

	case QualityCompressionParams::MET_BINARY:
	{
		ASSERT(binConfig.quaParams.binaryThreshold > 0 && binConfig.quaParams.binaryThreshold < 64);
		ASSERT(qbits == 1);

		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 c = quaReader_.GetBit();
			rec_.qua[i] = (char)(binConfig.archiveType.qualityOffset + (c ? QualityCompressionParams::Default::MaxThresholdValue
													: QualityCompressionParams::Default::MinThresholdValue));
		}
		break;
	}

	case QualityCompressionParams::MET_8BIN:
	{
		for (uint32 i = 0; i < rec_.seqLen; ++i)
		{
			uint32 c = quaReader_.GetBits(qbits);
			ASSERT(c < 8);
			rec_.qua[i] = (char)(binConfig.archiveType.qualityOffset + idxToQua_8bin[c]);
		}
		break;
	}

	case QualityCompressionParams::MET_QVZ:
	{
        for (uint32 i = 0; i < rec_.seqLen; ++i)
        {
            uint32 c = quaReader_.GetBits(qbits);
			rec_.qua[i] = (char)(c + binConfig.archiveType.qualityOffset);
        }
        break;
	}

	}
}

void IFastqPacker::ReadHeader(BitMemoryReader &/*metaReader_*/,
							  BitMemoryReader &headReader_,
							  const BinPackSettings &/*settings_*/,
							  FastqRecord &rec_)
{
	// now just store the header as it is -- TODO: more advanced compression
	ASSERT(rec_.head != NULL);

	rec_.headLen = headReader_.GetBits(8);			// max 255 length
	ASSERT(rec_.headLen > 0);

	// DEBUG
	std::fill(rec_.head, rec_.head + rec_.headLen, 'X');

	rec_.head[0] = '@';
	for (uint32 i = 1; i < rec_.headLen; ++i)
	{
		rec_.head[i] = headReader_.GetBits(7);		// max 128 ASCII values
	}
}



// TODO: divide and templatize
//
void FastqRecordsPackerSE::PackToBins(const std::map<uint32, FastqRecordsPtrBin>& dnaBins_,
									  BinaryBinBlock& binBlock_)
{
	binBlock_.Clear();
	binBlock_.blockType = BinaryBinBlock::MultiSignatureType;

	BitMemoryWriter metaWriter(binBlock_.metaData);
	BitMemoryWriter dnaWriter(binBlock_.dnaData);
	BitMemoryWriter quaWriter(binBlock_.quaData);
	BitMemoryWriter headWriter(binBlock_.headData);			// HINT: should be optional

	// store standard bins
	//
	const uint32 nBinId = binConfig.minimizer.TotalMinimizersCount();
	const auto nBinIter = dnaBins_.find(nBinId);

	uint64 totalRecordsCount = 0;
	for (auto i = dnaBins_.cbegin(); i != nBinIter; ++i)
	{
		uint32 binId = i->first;
		const FastqRecordsPtrBin& curBin = i->second;

		// skip empty bins
		ASSERT(curBin.records.size() != 0);
		ASSERT(binId != 0);		// DEBUG

		BinaryBinDescriptor& desc = binBlock_.descriptors[binId];
		desc.Clear();

#if DEV_DEBUG_MODE
		char minString[64];
		binConfig.minimizer.GenerateMinimizer(binId, minString);

		for (const auto r : i->second.records)
		{
			const char* minSeq = std::search(r->seq,
											 r->seq + r->seqLen,
											 minString,
											 minString + binConfig.minimizer.signatureLen);

			ASSERT(minSeq != r->seq + r->seqLen);
			ASSERT(minSeq == r->seq + r->minimPos);
		}
#endif

		PackToBin(curBin, metaWriter, dnaWriter, quaWriter, headWriter, desc, false);

		binBlock_.rawDnaSize += desc.rawDnaSize;
		binBlock_.rawHeadSize += desc.rawHeadSize;

		totalRecordsCount += desc.recordsCount;
	}

	// pack last N bin -- if present
	//
	if (nBinIter != dnaBins_.cend() && dnaBins_.at(nBinId).records.size() > 0)
	{
		BinaryBinDescriptor& nDesc = binBlock_.descriptors[nBinId];
		nDesc.Clear();

		const auto& nBin = dnaBins_.at(nBinId);

		PackToBin(nBin, metaWriter, dnaWriter, quaWriter, headWriter, nDesc, true);

		binBlock_.rawDnaSize += nDesc.rawDnaSize;
		binBlock_.rawHeadSize += nDesc.rawHeadSize;

		totalRecordsCount += nDesc.recordsCount;
	}

	binBlock_.metaSize = metaWriter.Position();
	binBlock_.dnaSize = dnaWriter.Position();
	binBlock_.quaSize = quaWriter.Position();
	binBlock_.headSize = headWriter.Position();
}




void FastqRecordsPackerSE::PackToBin(const std::vector<FastqRecord>& records_,
									 BinaryBinBlock& binBlock_,
									 uint32 binId_)
{
	ASSERT(binId_ == binConfig.minimizer.SignatureN());

	binBlock_.Clear();
	binBlock_.blockType = BinaryBinBlock::SingleSignatureType;

	BitMemoryWriter metaWriter(binBlock_.metaData);
	BitMemoryWriter dnaWriter(binBlock_.dnaData);
	BitMemoryWriter quaWriter(binBlock_.quaData);
	BitMemoryWriter headWriter(binBlock_.headData);		// HINT: should be optional


	// create temporary bin to store the data (as a workaround)
	// TODO: temporary solution, refactor
	FastqRecordsPtrBin curBin;
	for (const FastqRecord& rec : records_)
		curBin.records.push_back((FastqRecord*)&rec);

	curBin.stats.minSeqLen = curBin.stats.maxSeqLen= records_.front().seqLen;
	curBin.stats.minAuxLen = curBin.stats.maxAuxLen= records_.front().auxLen;

	binBlock_.auxDescriptors.resize(1);
	auto& curDesc = binBlock_.auxDescriptors.front();


	// store the data into binary form
	//
	PackToBin(curBin, metaWriter, dnaWriter, quaWriter, headWriter, curDesc, true);


	// update the block descriptor
	//
	binBlock_.rawDnaSize = curDesc.rawDnaSize;
	binBlock_.metaSize = curDesc.metaSize;
	binBlock_.dnaSize = curDesc.dnaSize;
	binBlock_.quaSize = curDesc.quaSize;
	binBlock_.signature = binId_;
	binBlock_.rawHeadSize = curDesc.rawHeadSize;
	binBlock_.headSize = curDesc.headSize;
}


void FastqRecordsPackerSE::PackToBin(const FastqRecordsPtrBin &fqBin_,
									BitMemoryWriter& metaWriter_,
									BitMemoryWriter& dnaWriter_,
									BitMemoryWriter& quaWriter_,
									BitMemoryWriter& headWriter_,
									BinaryBinDescriptor& binDesc_,
									bool nBin_)
{
	//const uint32 recordsCount = fqBin_.records.size();
	const uint64 initialMetaPos = metaWriter_.Position();
	const uint64 initialDnaPos = dnaWriter_.Position();
	const uint64 initialQuaPos = quaWriter_.Position();
	const uint64 initialHeadPos = headWriter_.Position();

	binDesc_.recordsCount = 0;

	ASSERT(fqBin_.records.size() != 0);
	//ASSERT(fqBin_.records.size() < (1 << 30));
	ASSERT(fqBin_.stats.minSeqLen > 0);
	ASSERT(fqBin_.stats.maxSeqLen > 0);
	ASSERT(fqBin_.stats.maxSeqLen >= fqBin_.stats.minSeqLen);

	BinPackSettings settings;
	settings.minLen = fqBin_.stats.minSeqLen;
	settings.maxLen = fqBin_.stats.maxSeqLen;
	settings.hasConstLen = (settings.minLen == settings.maxLen);
	settings.hasReadGroups = false;
	settings.usesHeaders = binConfig.archiveType.readsHaveHeaders;

	if (!nBin_)
		settings.suffixLen = binConfig.minimizer.signatureLen;
	else
		settings.suffixLen = 0;

	// if needed, the lenghts can be Huffman'd
	if (!settings.hasConstLen)
		settings.bitsPerLen = bit_length(fqBin_.stats.maxSeqLen - fqBin_.stats.minSeqLen);

	// store meta data header
	//
	metaWriter_.PutBits(fqBin_.stats.minSeqLen, LenBits);
	metaWriter_.PutBits(fqBin_.stats.maxSeqLen, LenBits);
	metaWriter_.PutBit(settings.hasReadGroups);


	// start packing records
	//
	StoreRecords(fqBin_.records, settings, metaWriter_, dnaWriter_, quaWriter_, headWriter_, binDesc_);


	// end packing
	//
	metaWriter_.FlushPartialWordBuffer();
	dnaWriter_.FlushPartialWordBuffer();
	quaWriter_.FlushPartialWordBuffer();
	headWriter_.FlushPartialWordBuffer();

	binDesc_.metaSize = metaWriter_.Position() - initialMetaPos;
	binDesc_.dnaSize = dnaWriter_.Position() - initialDnaPos;
	binDesc_.quaSize = quaWriter_.Position() - initialQuaPos;
	binDesc_.headSize = headWriter_.Position() - initialHeadPos;
}


void FastqRecordsPackerSE::UnpackFromBin(const BinaryBinBlock& binBin_,
									   std::vector<FastqRecord>& reads_,
									   FastqChunk& fqChunk_,
									   bool append_)
{
	ASSERT(binBin_.blockType == BinaryBinBlock::SingleSignatureType);
	ASSERT(binBin_.auxDescriptors.size() > 0);
	ASSERT(binBin_.signature != 0);


	// calculate global bin properties
	//
	uint64 recordsCount = 0;

	for (const auto& desc : binBin_.auxDescriptors)
	{
		ASSERT(desc.recordsCount > 0);
		recordsCount += desc.recordsCount;
	}

	//ASSERT(recordsCount < (1 << 30));
	ASSERT(recordsCount != 0);

	// initialize buffers
	//
	uint64 recIdx = 0;
	if (append_)
	{
		recIdx = reads_.size();
		recordsCount += recIdx;

		ASSERT(fqChunk_.data.Size() >= fqChunk_.size + binBin_.rawDnaSize * 2 + binBin_.rawHeadSize);
	}
	else
	{
		if (fqChunk_.data.Size() < binBin_.rawDnaSize * 2 + binBin_.rawHeadSize)
			fqChunk_.data.Extend(binBin_.rawDnaSize * 2 + binBin_.rawHeadSize);
	}

	// reset stats for initialization
	//
	if (recIdx == 0)
	{
		fqChunk_.size = 0;
		reads_.clear();

#if EXTRA_MEM_OPT
		reads_.shrink_to_fit();
#endif
	}

	reads_.resize(recordsCount);


	// initialize IO
	//
	BitMemoryReader metaReader(binBin_.metaData, binBin_.metaSize);
	BitMemoryReader dnaReader(binBin_.dnaData, binBin_.dnaSize);
	BitMemoryReader quaReader(binBin_.quaData, binBin_.quaSize);
	BitMemoryReader headReader(binBin_.headData, binBin_.headSize);


	// initialize settings
	//
	BinPackSettings settings;
	settings.signatureId = binBin_.signature;

	if (binBin_.signature != binConfig.minimizer.TotalMinimizersCount())	// N bin?
	{
		settings.suffixLen = binConfig.minimizer.signatureLen;
		binConfig.minimizer.GenerateMinimizer(binBin_.signature, settings.signatureString.data());
	}
	else
	{
		settings.suffixLen = 0;
	}
	settings.usesHeaders = binConfig.archiveType.readsHaveHeaders;


	// unpack records
	//
	for (const BinaryBinDescriptor& desc : binBin_.auxDescriptors)
	{
		ASSERT(desc.recordsCount > 0);

		const uint64 initialMetaPos = metaReader.Position();
		const uint64 initialDnaPos = dnaReader.Position();
		const uint64 initialQuaPos = quaReader.Position();
		const uint64 initialHeadPos = headReader.Position();

		// read bin header
		//
		settings.minLen = metaReader.GetBits(LenBits);
		settings.maxLen = metaReader.GetBits(LenBits);
		ASSERT(settings.minLen > 0);
		ASSERT(settings.maxLen >= settings.minLen);

		settings.hasReadGroups = metaReader.GetBit() != 0;
		ASSERT(!settings.hasReadGroups);

		settings.hasConstLen = (settings.minLen == settings.maxLen);
		ASSERT(settings.hasConstLen);

		if (!settings.hasConstLen)
		{
			settings.bitsPerLen = bit_length(settings.maxLen - settings.minLen);
			ASSERT(settings.bitsPerLen > 0);
		}

		// read the records contained in the block
		//
		ReadRecords(reads_.begin() + recIdx, reads_.begin() + recIdx + desc.recordsCount,
					settings, metaReader, dnaReader, quaReader, headReader, fqChunk_);

		recIdx += desc.recordsCount;

		metaReader.FlushInputWordBuffer();
		dnaReader.FlushInputWordBuffer();
		quaReader.FlushInputWordBuffer();
		headReader.FlushInputWordBuffer();

		ASSERT(metaReader.Position() - initialMetaPos == desc.metaSize);
		ASSERT(dnaReader.Position() - initialDnaPos == desc.dnaSize);
		ASSERT(quaReader.Position() - initialQuaPos == desc.quaSize);
		ASSERT(headReader.Position() - initialHeadPos == desc.headSize);
	}
}


void FastqRecordsPackerSE::StoreRecords(const std::vector<FastqRecord*>& records_,
										const BinPackSettings& settings_,
										BitMemoryWriter& metaWriter_,
										BitMemoryWriter& dnaWriter_,
										BitMemoryWriter& quaWriter_,
										BitMemoryWriter& headWriter_,
										BinaryBinDescriptor& binDesc_)
{
	for (auto rec : records_)
	{
		// handle the read length
		//
		if (!settings_.hasConstLen)
			metaWriter_.PutBits(rec->seqLen - settings_.minLen, settings_.bitsPerLen);

		// store single read
		//
		StoreNextRecord(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, *rec);
		binDesc_.rawDnaSize += rec->seqLen;
		binDesc_.rawHeadSize += rec->headLen;

		ASSERT(!binConfig.archiveType.readsHaveHeaders || rec->headLen > 0);

		binDesc_.recordsCount++;
	}
}


void FastqRecordsPackerSE::ReadRecords(std::vector<FastqRecord>::iterator itBegin_,
										std::vector<FastqRecord>::iterator itEnd_,
										const BinPackSettings& settings_,
										BitMemoryReader& metaReader_,
										BitMemoryReader& dnaReader_,
										BitMemoryReader& quaReader_,
										BitMemoryReader& headReader_,
										FastqChunk& fqChunk_)
{
	for ( ; itBegin_ != itEnd_; ++itBegin_)
	{
		FastqRecord& rec = *itBegin_;

		// reset the read and initialize the data pointers
		//
		rec.Reset();

		// get the read length
		//
		if (settings_.hasConstLen)
			rec.seqLen = settings_.minLen;
		else
			rec.seqLen = metaReader_.GetBits(settings_.bitsPerLen) + settings_.minLen;

		ASSERT(rec.seqLen > 0 && rec.seqLen < FastqRecord::MaxSeqLen);


		// bind sequence (and quality) for _1
		//
		rec.seq = (char*)(fqChunk_.data.Pointer() + fqChunk_.size);
		rec.qua = rec.seq + rec.seqLen;

		if (binConfig.archiveType.readsHaveHeaders)
		{
			rec.head = rec.qua + rec.seqLen;
		}


		bool b = ReadNextRecord(metaReader_, dnaReader_, quaReader_, headReader_, settings_, rec);
		ASSERT(b);

		if (settings_.suffixLen > 0)
		{
			std::copy(settings_.signatureString.data(),
					  settings_.signatureString.data() + binConfig.minimizer.signatureLen,
					  rec.seq + rec.minimPos);
		}

		fqChunk_.size += rec.seqLen * 2 + rec.headLen;		// dna + qua + head
	}
}


void FastqRecordsPackerPE::StoreRecords(const std::vector<FastqRecord*>& records_,
										const BinPackSettings& settings_,
										BitMemoryWriter& metaWriter_,
										BitMemoryWriter& dnaWriter_,
										BitMemoryWriter& quaWriter_,
										BitMemoryWriter& headWriter_,
										BinaryBinDescriptor& binDesc_)
{
	BinPackSettings pairSettings = settings_;
	pairSettings.suffixLen = 0;
	pairSettings.usesHeaders = false;

	for (const FastqRecord* rec : records_)
	{
		// store extra metadata information
		//
		if (!settings_.hasConstLen)
		{
			metaWriter_.PutBits(rec->seqLen - settings_.minLen, settings_.bitsPerLen);
			metaWriter_.PutBits(rec->auxLen - settings_.minLen, settings_.bitsPerLen);
		}

		if (settings_.suffixLen != 0)
			metaWriter_.PutBit(rec->IsPairSwapped());


		// store read _1
		//
		StoreNextRecord(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, *rec);
		binDesc_.rawDnaSize += rec->seqLen;


		// store read _2
		//
		FastqRecord pair = rec->GetPair();

		StoreNextRecord(metaWriter_, dnaWriter_, quaWriter_, headWriter_, pairSettings, pair);
		binDesc_.rawDnaSize += rec->auxLen;

		// WARN: only one copy of header
		binDesc_.rawHeadSize += rec->headLen;

		binDesc_.recordsCount++;
	}
}


void FastqRecordsPackerPE::ReadRecords(std::vector<FastqRecord>::iterator itBegin_,
										std::vector<FastqRecord>::iterator itEnd_,
										const BinPackSettings& settings_,
										BitMemoryReader& metaReader_,
										BitMemoryReader& dnaReader_,
										BitMemoryReader& quaReader_,
										BitMemoryReader& headReader_,
										FastqChunk& fqChunk_)
{
	BinPackSettings pairSettings = settings_;
	pairSettings.suffixLen = 0;
	pairSettings.usesHeaders = false;

	for ( ; itBegin_ != itEnd_; ++itBegin_)
	{
		FastqRecord& rec = *itBegin_;

		// reset the read and initialize the data pointers
		//
		rec.Reset();

		// read metadata information
		//
		if (settings_.hasConstLen)
		{
			rec.seqLen = settings_.minLen;
			rec.auxLen = settings_.minLen;
		}
		else
		{
			rec.seqLen = metaReader_.GetBits(settings_.bitsPerLen) + settings_.minLen;
			rec.auxLen = metaReader_.GetBits(settings_.bitsPerLen) + settings_.minLen;
		}

		if (settings_.suffixLen != 0)
			rec.SetPairSwapped(metaReader_.GetBit() != 0);

		// set-up read data pointers
		//
		rec.seq = (char*)(fqChunk_.data.Pointer() + fqChunk_.size);
		rec.qua = rec.seq + rec.seqLen + rec.auxLen;

		if (binConfig.archiveType.readsHaveHeaders)
		{
			rec.head = rec.qua + rec.seqLen + rec.auxLen;
		}


		// read _1
		//
		bool b = ReadNextRecord(metaReader_, dnaReader_, quaReader_, headReader_, settings_, rec);
		ASSERT(b);

		if (settings_.suffixLen > 0)
		{
			std::copy(settings_.signatureString.data(),
					  settings_.signatureString.data() + binConfig.minimizer.signatureLen,
					  rec.seq + rec.minimPos);
		}


		// read _2
		//
		FastqRecord pair = rec.GetPair();
		b = ReadNextRecord(metaReader_, dnaReader_, quaReader_, headReader_, pairSettings, pair);
		ASSERT(b);


		// update the chunk size
		//
		fqChunk_.size += (rec.seqLen + rec.auxLen) * 2 + rec.headLen;	// includes quality
	}
}
