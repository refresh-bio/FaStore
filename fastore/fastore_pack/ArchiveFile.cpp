/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../core/Globals.h"

#include <string>

#include "ArchiveFile.h"
#include "CompressedBlockData.h"

#include "../core/Exception.h"


ArchiveFileWriter::ArchiveFileWriter()
	:	metaStream(NULL)
	,	dataStream(NULL)
{}


ArchiveFileWriter::~ArchiveFileWriter()
{
	if (metaStream != NULL)
		delete metaStream;

	if (dataStream != NULL)
		delete dataStream;
}


void ArchiveFileWriter::StartCompress(const std::string &fileName_, const ArchiveConfig& config_)
{
	ASSERT(metaStream == NULL);
	ASSERT(dataStream == NULL);

	metaStream = new FileStreamWriter(fileName_ + ".cmeta");
	metaStream->SetBuffering(true);

	dataStream = new FileStreamWriter(fileName_ + ".cdata");


	// clear header and footer
	//
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(ArchiveFileHeader), 0);

	fileFooter.blockSizes.clear();
	fileFooter.signatures.clear();
	fileFooter.config = config_;


	// skip header pos
	//
	metaStream->SetPosition(ArchiveFileHeader::HeaderSize);
}


void ArchiveFileWriter::WriteNextBin(const DataChunk& compData_, uint32 signature_)
{
	ASSERT(compData_.size > 0);

	fileFooter.blockSizes.push_back(compData_.size);
	fileFooter.signatures.push_back(signature_);

	dataStream->Write(compData_.data.Pointer(), compData_.size);
}


void ArchiveFileWriter::FinishCompress()
{
	ASSERT(metaStream != NULL);
	ASSERT(dataStream != NULL);

	// prepare header and write footer
	//
	fileHeader.footerOffset = metaStream->Position();

	WriteFileFooter();


	// fill header and write
	//
	fileHeader.footerSize = metaStream->Position() - fileHeader.footerOffset;

	metaStream->SetPosition(0);
	WriteFileHeader();


	// cleanup exit
	//
	metaStream->Close();
	delete metaStream;
	metaStream = NULL;

	dataStream->Close();
	delete dataStream;
	dataStream = NULL;
}


void ArchiveFileWriter::WriteFileHeader()
{
	metaStream->Write((byte*)&fileHeader, ArchiveFileHeader::HeaderSize);
}


void ArchiveFileWriter::WriteFileFooter()
{
	// store the archive blocks info
	//
	uint32 blockCount = fileFooter.blockSizes.size();
	metaStream->Write((byte*)&blockCount, sizeof(uint32));
	metaStream->Write((byte*)fileFooter.blockSizes.data(), fileFooter.blockSizes.size() * sizeof(uint64));
	metaStream->Write((byte*)fileFooter.signatures.data(), fileFooter.signatures.size() * sizeof(uint32));


	// store the archive configuration
	//
	metaStream->Write((byte*)&fileFooter.config, sizeof(ArchiveConfig));


	if (fileFooter.config.quaParams.method != QualityCompressionParams::MET_QVZ
			&& !fileFooter.config.archType.readsHaveHeaders)
		return;


	// store the quality compression data
	//
	Buffer mem(2 << 10);
	BitMemoryWriter writer(mem);

	if (fileFooter.config.quaParams.method == QualityCompressionParams::MET_QVZ)
	{
		// store well rng
		writer.PutBytes((byte*)fileFooter.quaData.well.state, sizeof(fileFooter.quaData.well.state));

		// store max read length -- the max number of columns
		writer.PutBytes((byte*)&fileFooter.quaData.max_read_length, sizeof(fileFooter.quaData.max_read_length));

		// store codebook
		fileFooter.quaData.codebook.WriteCodebook(writer, fileFooter.quaData.max_read_length);
	}


	// store the read headers compression stuff -- based on BinFile, TODO: one serialization method
	//
	if (fileFooter.config.archType.readsHaveHeaders)
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

				ASSERT(f.possibleValues.size() < 255);			// we shall set max limit over the number of available tokens

				if (!f.isConst)
				{
					writer.PutByte(f.possibleValues.size());
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
		if (fileFooter.config.archType.readType == ArchiveType::READ_PE)
		{
			// HINT: may also lack the identifier
			//ASSERT(fileFooter.headData.pairedEndFieldIdx != 0);

			writer.PutByte(fileFooter.headData.pairedEndFieldIdx);
		}
	}


	writer.FlushPartialWordBuffer();
	metaStream->Write(writer.Pointer(), writer.Position());

}


ArchiveFileReader::ArchiveFileReader()
	:	metaStream(NULL)
	,	dataStream(NULL)
	//,	blockIdx(0)
{}


ArchiveFileReader::~ArchiveFileReader()
{
	if (metaStream != NULL)
		delete metaStream;

	if (dataStream != NULL)
		delete dataStream;
}


void ArchiveFileReader::StartDecompress(const std::string &fileName_, ArchiveConfig& config_)
{
	ASSERT(metaStream == NULL);

	metaStream = new FileStreamReader(fileName_ + ".cmeta");
	metaStream->SetBuffering(true);

	dataStream = new FileStreamReader(fileName_ + ".cdata");

	if (metaStream->Size() == 0 || dataStream->Size() == 0)
		throw Exception("Empty archive.");

	// Read file header
	//
	std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(ArchiveFileHeader), 0);
	ReadFileHeader();

	if (fileHeader.footerOffset + (uint64)fileHeader.footerSize > metaStream->Size())
	{
		delete metaStream;
		metaStream = NULL;
		throw Exception("Corrupted archive.");
	}

	// clean footer
	//
	fileFooter.blockSizes.clear();

	metaStream->SetPosition(fileHeader.footerOffset);
	ReadFileFooter();

	metaStream->SetPosition(ArchiveFileHeader::HeaderSize);
	config_ = fileFooter.config;


	// initialize block iterator
	//
	uint64 offset = 0;
	blockArray.resize(fileFooter.blockSizes.size());
	for (uint64 i = 0; i < fileFooter.blockSizes.size(); ++i)
	{
		blockArray[i] = std::make_tuple(fileFooter.signatures[i], fileFooter.blockSizes[i], offset);
		offset +=  fileFooter.blockSizes[i];
	}
//	std::random_shuffle(blockArray.begin(), blockArray.end());
	blockIterator = blockArray.begin();
}


void ArchiveFileReader::ReadFileHeader()
{
	metaStream->Read((byte*)&fileHeader, ArchiveFileHeader::HeaderSize);
}


void ArchiveFileReader::ReadFileFooter()
{
	// read the archive blocks information
	//
	uint32 blockCount = 0;
	metaStream->Read((byte*)&blockCount, sizeof(uint32));
	ASSERT(blockCount > 0);
	fileFooter.blockSizes.resize(blockCount);
	fileFooter.signatures.resize(blockCount);

	metaStream->Read((byte*)fileFooter.blockSizes.data(), fileFooter.blockSizes.size() * sizeof(uint64));
	metaStream->Read((byte*)fileFooter.signatures.data(), fileFooter.signatures.size() * sizeof(uint32));


	// read the archive config
	//
	metaStream->Read((byte*)&fileFooter.config, sizeof(ArchiveConfig));


	// read quality compression data
	//
	const uint64 extraDataSize = metaStream->Size() - metaStream->Position();
	if (extraDataSize == 0)
		return;

	Buffer mem(extraDataSize);
	metaStream->Read(mem.Pointer(), extraDataSize);
	BitMemoryReader reader(mem, extraDataSize);

	if (fileFooter.config.quaParams.method == QualityCompressionParams::MET_QVZ)
	{
		// read well
		reader.GetBytes((byte*)fileFooter.quaData.well.state, sizeof(fileFooter.quaData.well.state));

		// read the max read length -- the max number of columns for qvz
		reader.GetBytes((byte*)&fileFooter.quaData.max_read_length, sizeof(fileFooter.quaData.max_read_length));
		ASSERT(fileFooter.quaData.max_read_length > 0 && fileFooter.quaData.max_read_length < FastqRecord::MaxSeqLen);

		// read codebook
        struct alphabet_t *A = alloc_alphabet(ALPHABET_SIZE);
		fileFooter.quaData.codebook.ReadCodebook(reader, A, fileFooter.quaData.max_read_length);
	}


	// read headers stuff -- copied from BinFile
	//
	if (fileFooter.config.archType.readsHaveHeaders)
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
					possibleValues = reader.GetByte();
					ASSERT(possibleValues > 1 && possibleValues < 255);			// we shall set max limit over the number of available tokens
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
		if (fileFooter.config.archType.readType == ArchiveType::READ_PE)
		{
			fileFooter.headData.pairedEndFieldIdx = reader.GetByte();

			// HINT: may also lack this identifier
			//ASSERT(fileFooter.headData.pairedEndFieldIdx != 0);

			ASSERT(fileFooter.headData.pairedEndFieldIdx < fileFooter.headData.fields.size());
		}
	}
}


bool ArchiveFileReader::ReadNextBin(DataChunk& buffer_, uint32& signature_)
{
#if 0
	if (blockIdx >= fileFooter.blockSizes.size())
	{
		signature_ = 0;
		buffer_.size = 0;
		return false;
	}

	const uint64 bs = fileFooter.blockSizes[blockIdx];
	if (buffer_.data.Size() < bs)
		buffer_.data.Extend(bs + (bs / 8));

	dataStream->Read(buffer_.data.Pointer(), bs);
	buffer_.size = bs;
	signature_ = fileFooter.signatures[blockIdx];

	blockIdx++;
#else

	if (blockIterator == blockArray.end())
	{
		signature_ = 0;
		buffer_.size = 0;
		return false;
	}

	const uint64 bs = std::get<1>(*blockIterator);
	if (buffer_.data.Size() < bs)
		buffer_.data.Extend(bs + (bs / 8));

	dataStream->SetPosition(std::get<2>(*blockIterator));

	dataStream->Read(buffer_.data.Pointer(), bs);
	buffer_.size = bs;
	signature_ = std::get<0>(*blockIterator);

	blockIterator++;

#endif
	return true;
}


void ArchiveFileReader::FinishDecompress()
{
	ASSERT(metaStream != NULL);
	ASSERT(dataStream != NULL);

	metaStream->Close();
	delete metaStream;
	metaStream = NULL;

	dataStream->Close();
	delete dataStream;
	dataStream = NULL;
}
