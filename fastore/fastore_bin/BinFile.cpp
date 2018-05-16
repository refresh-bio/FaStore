/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"

#include <cstring>

#include "BinFile.h"
#include "BitMemory.h"
#include "BinBlockData.h"
#include "Exception.h"
#include "BitMemory.h"
#include "Utils.h"
#include "QVZ.h"


BinFileWriter::BinFileWriter()
	:	metaStream(NULL)
	,	dnaStream(NULL)
	,	quaStream(NULL)
	,	headStream(NULL)
{}


BinFileWriter::~BinFileWriter()
{
	if (metaStream != NULL)
		delete metaStream;

	if (dnaStream != NULL)
		delete dnaStream;

	if (quaStream != NULL)
		delete quaStream;

	if (headStream != NULL)
		delete headStream;
}


void BinFileWriter::StartCompress(const std::string& fileName_, const BinModuleConfig& params_)
{
	ASSERT(metaStream == NULL);
	ASSERT(dnaStream == NULL);

	metaStream = new FileStreamWriter(fileName_ + ".bmeta");
	((FileStreamWriter*)metaStream)->SetBuffering(true);

	dnaStream = new FileStreamWriter(fileName_ + ".bdna");
	((FileStreamWriter*)dnaStream)->SetBuffering(true);

	quaStream = new FileStreamWriter(fileName_ + ".bqua");
	((FileStreamWriter*)dnaStream)->SetBuffering(true);

	if (params_.archiveType.readsHaveHeaders)
	{
		headStream = new FileStreamWriter(fileName_ + ".bhead");
		((FileStreamWriter*)headStream)->SetBuffering(true);
	}


	// clear header and footer
	//
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(BinFileHeader), 0);

	fileFooter.Clear();
	fileFooter.params = params_;

	// skip header pos
	//
	metaStream->SetPosition(BinFileHeader::HeaderSize);

	// clear stats
	//
	globalFastqStats.Clear();
}


void BinFileWriter::WriteNextBlock(const BinaryBinBlock* block_)
{
	ASSERT(block_ != NULL);

	// TODO: better create an inherited virtual class to handle this case!!!!
	if (block_->blockType == BinaryBinBlock::MultiSignatureType)
	{
		ASSERT(block_->descriptors.size() > 0);

		uint64 metaOffset = 0;
		uint64 dnaOffset = 0;
		uint64 quaOffset = 0;
		uint64 headOffset = 0;

		for (const auto& iDesc : block_->descriptors)
		{
			// update stats
			//
			fileHeader.recordsCount += iDesc.second.recordsCount;

			BinFileFooter::BinInfo& fo = fileFooter.binOffsets[iDesc.first];

			BinFileFooter::BlockMetaData bmd(iDesc.second);
			bmd.metaFileOffset = metaStream->Position();
			bmd.dnaFileOffset = dnaStream->Position();
			bmd.quaFileOffset = quaStream->Position();

			fo.totalMetaSize += bmd.metaSize;
			fo.totalDnaSize += bmd.dnaSize;
			fo.totalQuaSize += bmd.quaSize;

			fo.totalRawDnaSize += bmd.rawDnaSize;
			fo.totalRecordsCount += bmd.recordsCount;


			// store block
			//
			metaStream->Write(block_->metaData.Pointer() + metaOffset, bmd.metaSize);
			dnaStream->Write(block_->dnaData.Pointer() + dnaOffset, bmd.dnaSize);
			quaStream->Write(block_->quaData.Pointer() + quaOffset, bmd.quaSize);

			metaOffset += bmd.metaSize;
			dnaOffset += bmd.dnaSize;
			quaOffset += bmd.quaSize;


			if (fileFooter.params.archiveType.readsHaveHeaders)
			{
				ASSERT(bmd.headSize > 0);
				ASSERT(bmd.rawHeadSize > 0);

				bmd.headFileOffset = headStream->Position();

				fo.totalRawHeadSize += bmd.rawHeadSize;
				fo.totalHeadSize += bmd.headSize;

				headStream->Write(block_->headData.Pointer() + headOffset, bmd.headSize);

				headOffset += bmd.headSize;
			}

			fo.blocksMetaData.push_back(bmd);
		}

		ASSERT(metaOffset == block_->metaSize);
		ASSERT(dnaOffset == block_->dnaSize);
		ASSERT(quaOffset == block_->quaSize);
		ASSERT(headOffset == block_->headSize);
	}
	else
	{
		ASSERT(block_->auxDescriptors.size() > 0);
		ASSERT(block_->signature != 0);

		BinFileFooter::BinInfo& fo = fileFooter.binOffsets[block_->signature];

		uint64 metaOffset = 0;
		uint64 dnaOffset = 0;
		uint64 quaOffset = 0;
		uint64 headOffset = 0;
		for (const auto& iDesc : block_->auxDescriptors)
		{
			// update stats
			//
			fileHeader.recordsCount += iDesc.recordsCount;

			BinFileFooter::BlockMetaData bmd(iDesc);
			bmd.metaFileOffset = metaStream->Position();
			bmd.dnaFileOffset = dnaStream->Position();
			bmd.quaFileOffset = quaStream->Position();

			fo.totalMetaSize += bmd.metaSize;
			fo.totalDnaSize += bmd.dnaSize;
			fo.totalQuaSize += bmd.quaSize;

			fo.totalRawDnaSize += bmd.rawDnaSize;
			fo.totalRecordsCount += bmd.recordsCount;

			// store block
			//
			metaStream->Write(block_->metaData.Pointer() + metaOffset, bmd.metaSize);
			dnaStream->Write(block_->dnaData.Pointer() + dnaOffset, bmd.dnaSize);
			quaStream->Write(block_->quaData.Pointer() + quaOffset, bmd.quaSize);

			metaOffset += bmd.metaSize;
			dnaOffset += bmd.dnaSize;
			quaOffset += bmd.quaSize;

			if (fileFooter.params.archiveType.readsHaveHeaders)
			{
				ASSERT(bmd.headSize > 0);
				ASSERT(bmd.rawHeadSize > 0);

				bmd.headFileOffset = headStream->Position();

				fo.totalRawHeadSize += bmd.rawHeadSize;
				fo.totalHeadSize += bmd.headSize;

				headStream->Write(block_->headData.Pointer() + headOffset, bmd.headSize);

				headOffset += bmd.headSize;

			}

			fo.blocksMetaData.push_back(bmd);
		}

		ASSERT(metaOffset == block_->metaSize);
		ASSERT(dnaOffset == block_->dnaSize);
		ASSERT(quaOffset == block_->quaSize);
		ASSERT(headOffset == block_->headSize);
	}


	// update stats
	//
	globalFastqStats.Update(block_->stats);
}


void BinFileWriter::FinishCompress()
{
	ASSERT(metaStream != NULL);
	ASSERT(dnaStream != NULL);
	ASSERT(quaStream != NULL);


	// prepare header and footer
	//
	std::fill(fileHeader.reserved, fileHeader.reserved + BinFileHeader::ReservedBytes, 0);
	fileHeader.blockCount = fileFooter.binOffsets.size();
	fileHeader.footerOffset = metaStream->Position();
	fileHeader.usesHeaderStream = fileFooter.params.archiveType.readsHaveHeaders;

	// prepare quality stuff
	//
	if (fileFooter.params.quaParams.method == QualityCompressionParams::MET_QVZ
			&& fileFooter.params.binningLevel == 0)
	{
		// Compute the marginal pmfs of the global stats
		globalFastqStats.Compute_marginal_pmf();

		// Compute the QVZ quantizers
		fileFooter.quaData.codebook.ComputeFromStats(globalFastqStats.qua.training_stats,
													&fileFooter.params.quaParams.qvzOpts);

		// We need to initialize the WELL RND
		memset(&fileFooter.quaData.well, 0, sizeof(struct well_state_t));

		// Initialize WELL state vector with libc rand
		srand((uint32_t) time(0));
		for (uint32_t i = 0; i < 32; ++i)
		{
#ifndef DEBUG
		fileFooter.quaData.well.state[i] = rand();
#else
		fileFooter.quaData.well.state[i] = 0x55555555;
#endif
		}

		// set the maximum calculated read length for future allocation
		// of qvz columnar statistics
		fileFooter.quaData.max_read_length = MAX(globalFastqStats.maxSeqLen, globalFastqStats.maxAuxLen);
	}


	// prepare quality stuff
	//
	if (fileHeader.usesHeaderStream
			&& fileFooter.params.binningLevel == 0)
	{
		fileFooter.headData = globalFastqStats.head;
	}


	WriteFileFooter();
	fileHeader.footerSize = metaStream->Position() - fileHeader.footerOffset;

	metaStream->SetPosition(0);
	WriteFileHeader();


	// cleanup
	//
	metaStream->Close();
	dnaStream->Close();
	quaStream->Close();

	delete metaStream;
	delete dnaStream;
	delete quaStream;

	metaStream = NULL;
	dnaStream = NULL;
	quaStream = NULL;

	if (headStream != NULL)
	{
		headStream->Close();
		delete headStream;
		headStream = NULL;
	}
}


void BinFileWriter::WriteFileHeader()
{
	metaStream->Write((byte*)&fileHeader, BinFileHeader::HeaderSize);
}


void BinFileWriter::WriteFileFooter()
{
	const uint64 osize = BinFileFooter::ParametersSize
			+ fileFooter.binOffsets.size() * sizeof(uint64) * 2;

	Buffer obuf(osize);
	BitMemoryWriter writer(obuf);
	writer.PutBytes((byte*)&fileFooter.params, BinFileFooter::ParametersSize);

	// calculate bin occupancy bitmap
	//
	std::vector<bool> signatureBitmap;
	signatureBitmap.resize(fileFooter.params.minimizer.TotalMinimizersCount() + 1, false);

	for (auto off : fileFooter.binOffsets)
		signatureBitmap[off.first] = true;


	// write the signature bitmap
	// TODO: direct write
	//
	for (std::vector<bool>::const_iterator iSig = signatureBitmap.begin();
		 iSig != signatureBitmap.end();
		 iSig++)
	{
		writer.PutBit(*iSig);
	}
	writer.FlushPartialWordBuffer();


	// calculate sub-footer offsets and save
	//



	// write the file offsets
	// TODO: delta encode
	// TODO: direct file write of offsets
	//
	for (const auto& iOff : fileFooter.binOffsets)
	{
		// those folks can be now calculated later - we don't need to store them
		//
		writer.Put8Bytes(iOff.second.totalMetaSize);
		writer.Put8Bytes(iOff.second.totalDnaSize);
		writer.Put8Bytes(iOff.second.totalQuaSize);

		writer.Put8Bytes(iOff.second.totalRawDnaSize);
		writer.Put8Bytes(iOff.second.totalRecordsCount);
		//


		// headers stuff
		//
		if (fileHeader.usesHeaderStream)
		{
			writer.Put8Bytes(iOff.second.totalHeadSize);
			writer.Put8Bytes(iOff.second.totalRawHeadSize);
		}


		// TODO: compress the offsets
		//
		// TODO: do not store header offsets if not used
		//
		writer.Put8Bytes(iOff.second.blocksMetaData.size());
		writer.PutBytes((byte*)iOff.second.blocksMetaData.data(),
						iOff.second.blocksMetaData.size() * sizeof(BinFileFooter::BlockMetaData));
	}


	// store the quality compression stuff
	//
	if (fileFooter.params.quaParams.method == QualityCompressionParams::MET_QVZ)
	{		
		writer.PutBytes((byte*)fileFooter.quaData.well.state, sizeof(fileFooter.quaData.well.state));
		writer.PutBytes((byte*)&fileFooter.quaData.max_read_length, sizeof(fileFooter.quaData.max_read_length));
		fileFooter.quaData.codebook.WriteCodebook(writer, fileFooter.quaData.max_read_length);
	}


	// store the read headers compression stuff
	//
	if (fileHeader.usesHeaderStream)
	{
		ASSERT(fileFooter.headData.fields.size() > 0);
		writer.PutByte(fileFooter.headData.fields.size());

		for (const FastqRawBlockStats::HeaderStats::Field& f : fileFooter.headData.fields)
		{
			writer.PutByte(f.isNumeric);
			writer.PutByte(f.isConst);
			writer.PutByte(f.separator);
			if (f.isNumeric)
			{
				writer.Put8Bytes(f.minValue);
				if (!f.isConst)
				{
					writer.Put8Bytes(f.maxValue);
				}
			}
			else
			{
				// not storing min/max length

				ASSERT(f.possibleValues.size() < FastqRawBlockStats::HeaderStats::MaxPossibleValues);

				if (!f.isConst)
				{
					writer.Put2Bytes(f.possibleValues.size());
				}

				for (const std::string& s : f.possibleValues)
				{
					writer.PutByte(s.size());
					writer.PutBytes((byte*)s.c_str(), s.size());
				}
			}
		}


		// extra in PE mode -- check which field is PE identifier
		//
		if (fileFooter.params.archiveType.readType == ArchiveType::READ_PE)
		{
			fileFooter.headData.pairedEndFieldIdx = 0;
			for (int32 i = fileFooter.headData.fields.size() - 1; i >= 0; i--)
			{
				const auto& f = fileFooter.headData.fields[i];
				if (f.isNumeric)
				{
					if (f.minValue == 1 && f.maxValue == 2)
					{
						fileFooter.headData.pairedEndFieldIdx = i;
						break;
					}
				}
			}
			// HINT: may also lack the identifier
			//ASSERT(fileFooter.headData.pairedEndFieldIdx != 0);

			writer.PutByte(fileFooter.headData.pairedEndFieldIdx);
		}
	}

	metaStream->Write(writer.Pointer(), writer.Position());
}


BinFileReader::BinFileReader()
	:	metaStream(NULL)
	,	dnaStream(NULL)
	,	quaStream(NULL)
	,	headStream(NULL)
{
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(BinFileHeader), 0);
}


BinFileReader::~BinFileReader()
{
	if (metaStream)
		delete metaStream;

	if (dnaStream)
		delete dnaStream;

	if (quaStream)
		delete quaStream;

	if (headStream)
		delete headStream;
}


void BinFileReader::StartDecompress(const std::string& fileName_, BinModuleConfig& params_)
{
	ASSERT(metaStream == NULL);
	ASSERT(dnaStream == NULL);

	metaStream = new FileStreamReader(fileName_ + ".bmeta");
	((FileStreamReader*)metaStream)->SetBuffering(true);

	if (metaStream->Size() == 0)
		throw Exception("Empty file.");

	dnaStream = new FileStreamReader(fileName_ + ".bdna");
	//((FileStreamReader*)dnaStream)->SetBuffering(true);
	quaStream = new FileStreamReader(fileName_ + ".bqua");
	//((FileStreamReader*)quaStream)->SetBuffering(true);

	// read header
	//
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(BinFileHeader), 0);
	ReadFileHeader();

	if ((fileHeader.blockCount == 0ULL) || fileHeader.footerOffset + (uint64)fileHeader.footerSize > metaStream->Size())
	{
		delete metaStream;
		delete dnaStream;
		delete quaStream;
		metaStream = NULL;
		dnaStream = NULL;
		quaStream = NULL;
		throw Exception("Corrupted archive header");
	}

	if (fileHeader.usesHeaderStream)
	{
		headStream = new FileStreamReader(fileName_ + ".bhead");
	//	((FileStreamReader*)headStream )->SetBuffering(true);
	}


	// read footer
	//
	fileFooter.Clear();

	metaStream->SetPosition(fileHeader.footerOffset);
	ReadFileFooter();
	params_ = fileFooter.params;

	metaStream->SetPosition(BinFileHeader::HeaderSize);

	metaOffsetIterator = fileFooter.binOffsets.begin();
}


bool BinFileReader::ReadNextBlock(BinaryBinBlock* block_)
{
	if (metaOffsetIterator == fileFooter.binOffsets.end())
		return false;

	ReadBlock(metaOffsetIterator->first, block_);

	metaOffsetIterator++;
	return true;
}


void BinFileReader::ReadBlock(uint32 signature_, BinaryBinBlock* block_)
{
	ASSERT(block_ != NULL);
	ASSERT(fileFooter.binOffsets.count(signature_) != 0);

	block_->Clear();
	block_->blockType = BinaryBinBlock::SingleSignatureType;
	block_->signature = signature_;

	const BinFileFooter::BinInfo& footOff = fileFooter.binOffsets.at(signature_);
	ASSERT(footOff.totalMetaSize > 0);
	ASSERT(footOff.totalDnaSize > 0);
	ASSERT(footOff.totalQuaSize > 0);
	ASSERT(footOff.totalRawDnaSize > 0);
	ASSERT(footOff.totalRecordsCount > 0);

	if (block_->metaData.Size() < footOff.totalMetaSize)
		block_->metaData.Extend(footOff.totalMetaSize);

	if (block_->dnaData.Size() < footOff.totalDnaSize)
		block_->dnaData.Extend(footOff.totalDnaSize);

	if (block_->quaData.Size() < footOff.totalQuaSize)
		block_->quaData.Extend(footOff.totalQuaSize);

	if (fileHeader.usesHeaderStream)
	{
		if (block_->headData.Size() < footOff.totalHeadSize)
			block_->headData.Extend(footOff.totalHeadSize);
	}


	uint64 metaOffset = 0;
	uint64 dnaOffset = 0;
	uint64 quaOffset = 0;
	uint64 headOffset = 0;
	for (std::vector<BinFileFooter::BlockMetaData>::const_iterator iBlock = footOff.blocksMetaData.begin();
		 iBlock != footOff.blocksMetaData.end();
		 iBlock++)
	{
		metaStream->SetPosition(iBlock->metaFileOffset);
		dnaStream->SetPosition(iBlock->dnaFileOffset);
		quaStream->SetPosition(iBlock->quaFileOffset);

		metaStream->Read(block_->metaData.Pointer() + metaOffset, iBlock->metaSize);
		dnaStream->Read(block_->dnaData.Pointer() + dnaOffset, iBlock->dnaSize);
		quaStream->Read(block_->quaData.Pointer() + quaOffset, iBlock->quaSize);

		metaOffset += iBlock->metaSize;
		dnaOffset += iBlock->dnaSize;
		quaOffset += iBlock->quaSize;

		block_->auxDescriptors.push_back(*iBlock);

		block_->metaSize += iBlock->metaSize;
		block_->dnaSize += iBlock->dnaSize;
		block_->quaSize += iBlock->quaSize;
		block_->rawDnaSize += iBlock->rawDnaSize;

		if (fileHeader.usesHeaderStream)
		{
			headStream->SetPosition(iBlock->headFileOffset);
			headStream->Read(block_->headData.Pointer() + headOffset, iBlock->headSize);

			headOffset += iBlock->headSize;

			ASSERT(iBlock->headSize > 0);
			ASSERT(iBlock->rawHeadSize > 0);

			block_->headSize += iBlock->headSize;
			block_->rawHeadSize += iBlock->rawHeadSize;
		}
	}

	// now this signature associated data can be freed here
	//
}


void BinFileReader::FinishDecompress()
{
	if (metaStream)
	{
		metaStream->Close();
		delete metaStream;
		metaStream = NULL;
	}

	if (dnaStream)
	{
		dnaStream->Close();
		delete dnaStream;
		dnaStream = NULL;
	}

	if (quaStream)
	{
		quaStream->Close();
		delete quaStream;
		quaStream = NULL;
	}

	if (headStream)
	{
		headStream->Close();
		delete headStream;
		headStream = NULL;
	}
}


void BinFileReader::ReadFileHeader()
{
	metaStream->Read((byte*)&fileHeader, BinFileHeader::HeaderSize);
}


void BinFileReader::ReadFileFooter()
{
	Buffer buffer(fileHeader.footerSize);

	// TODO: just read partial data --> the file offsets read directly from
	// disk to avoid unnecessary copy
	//
	metaStream->Read(buffer.Pointer(), fileHeader.footerSize);

	BitMemoryReader reader(buffer, fileHeader.footerSize);

	reader.GetBytes((byte*)&fileFooter.params, BinFileFooter::ParametersSize);


	// calculate bin occupancy bitmap
	//
	std::vector<bool> signatureBitmap;
	signatureBitmap.resize(fileFooter.params.minimizer.TotalMinimizersCount() + 1, false);

	for (uint32 i = 0; i < fileFooter.params.minimizer.TotalMinimizersCount() + 1; ++i)
		signatureBitmap[i] = reader.GetBit() != 0;

	reader.FlushInputWordBuffer();


	// read file offsets
	//
	for (uint32 i = 0; i < fileFooter.params.minimizer.TotalMinimizersCount() + 1; ++i)
	{
		if (!signatureBitmap[i])
			continue;

		BinFileFooter::BinInfo& fo = fileFooter.binOffsets[i];

		fo.totalMetaSize = reader.Get8Bytes();
		ASSERT(fo.totalMetaSize  > 0);

		fo.totalDnaSize = reader.Get8Bytes();
		ASSERT(fo.totalDnaSize > 0);

		fo.totalQuaSize = reader.Get8Bytes();
		ASSERT(fo.totalQuaSize > 0);

		fo.totalRawDnaSize = reader.Get8Bytes();
		ASSERT(fo.totalRawDnaSize > 0);

		fo.totalRecordsCount = reader.Get8Bytes();
		ASSERT(fo.totalRecordsCount > 0);

		if (fileHeader.usesHeaderStream)
		{
			fo.totalHeadSize = reader.Get8Bytes();
			ASSERT(fo.totalHeadSize > 0);

			fo.totalRawHeadSize = reader.Get8Bytes();
			ASSERT(fo.totalRawHeadSize > 0);
		}

		uint64 fcc = reader.Get8Bytes();
		ASSERT(fcc > 0);
		fo.blocksMetaData.resize(fcc);
		reader.GetBytes((byte*)fo.blocksMetaData.data(),
						fcc * sizeof(BinFileFooter::BlockMetaData));
	}


	// read the quality compression stuff
	//
	if (fileFooter.params.quaParams.method == QualityCompressionParams::MET_QVZ)
	{
		reader.GetBytes((byte*)fileFooter.quaData.well.state, sizeof(fileFooter.quaData.well.state));

		reader.GetBytes((byte*)&fileFooter.quaData.max_read_length, sizeof(fileFooter.quaData.max_read_length));
		ASSERT(fileFooter.quaData.max_read_length > 0 && fileFooter.quaData.max_read_length < FastqRecord::MaxSeqLen);

		struct alphabet_t *A = alloc_alphabet(ALPHABET_SIZE);
        //TODO fix this. Should we write the read length in the footer also?
		fileFooter.quaData.codebook.ReadCodebook(reader, A, fileFooter.quaData.max_read_length);


		// shall we dealloc A???
        //not here... We needed it to write the codebook back into the file footer...
	}


	// read headers stuff -- TODO: move to external serialization fcn
	//
	if (fileHeader.usesHeaderStream)
	{
		uint32 fieldsCount = reader.GetByte();
		ASSERT(fieldsCount > 0);

		fileFooter.headData.fields.resize(fieldsCount);

		for (FastqRawBlockStats::HeaderStats::Field& f : fileFooter.headData.fields)
		{
			f.isNumeric = reader.GetByte() != 0;
			f.isConst = reader.GetByte();
			f.separator = reader.GetByte();
			if (f.isNumeric)
			{
				f.minValue = reader.Get8Bytes();
				if (!f.isConst)
				{
					f.maxValue = reader.Get8Bytes();
				}
			}
			else
			{
				// not reading min/max length

				uint32 possibleValues = 1;
				if (!f.isConst)
				{
					possibleValues = reader.Get2Bytes();
					ASSERT(possibleValues > 1 && possibleValues < FastqRawBlockStats::HeaderStats::MaxPossibleValues);
				}

				for (uint32 i = 0; i < possibleValues; ++i)
				{
					uint32 ss = reader.GetByte();
					char sbuf[ss];
					reader.GetBytes((byte*)sbuf, ss);
					f.possibleValues.insert(std::string(sbuf, ss));
				}
			}
		}

		// read info about the field number defining SE/PE
		//
		if (fileFooter.params.archiveType.readType == ArchiveType::READ_PE)
		{
			fileFooter.headData.pairedEndFieldIdx = reader.GetByte();

			// HINT: may also lack this identifier
			//ASSERT(fileFooter.headData.pairedEndFieldIdx != 0);

			ASSERT(fileFooter.headData.pairedEndFieldIdx < fileFooter.headData.fields.size());
		}
	}
}
