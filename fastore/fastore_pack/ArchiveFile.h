/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DNARCHFILE
#define H_DNARCHFILE

#include "../fastore_bin/Globals.h"

#include "Params.h"
#include "CompressedBlockData.h"

#include "../fastore_bin/FileStream.h"
#include "../fastore_bin/Params.h"

#include "../fastore_bin/QVZ.h"


/**
 * A compressed FaStore archive file
 *
 */
class IArchiveFile
{
public:
	virtual ~IArchiveFile() {}

	struct ArchiveConfig
	{
		ArchiveType archType;
		MinimizerParameters minParams;
		QualityCompressionParams quaParams;
	};

protected:
	struct ArchiveFileHeader
	{
		static const uint64 ReservedBytes = 8;
		static const uint64 HeaderSize = 8 + 8 + ReservedBytes;

		uint64 footerOffset;
		uint64 footerSize;

		uchar reserved[ReservedBytes];

		ArchiveFileHeader()
		{
			STATIC_ASSERT(sizeof(ArchiveFileHeader) == HeaderSize);
		}
	};

	struct ArchiveFileFooter
	{
		std::vector<uint64> blockSizes;		// TODO: compress space
		std::vector<uint32> signatures;		// TODO: this can be reduced to bitmap

		ArchiveConfig config;

		QualityCompressionData quaData;
		FastqRawBlockStats::HeaderStats headData;
	};

	ArchiveFileHeader fileHeader;
	ArchiveFileFooter fileFooter;
};


class ArchiveFileWriter : public IArchiveFile
{
public:
	ArchiveFileWriter();
	~ArchiveFileWriter();

	void StartCompress(const std::string& fileName_, const ArchiveConfig& config_);
    void WriteNextBin(const DataChunk& compData_, uint32 signature_);
	void FinishCompress();


	// TODO: update accordingly to QVZ required data for decompression
	//
	void SetQualityCompressionData(const QualityCompressionData& quaData_)
	{
		fileFooter.quaData = quaData_;
	}

	void SetHeadersCompressionData(const FastqRawBlockStats::HeaderStats& headData_)
	{
		fileFooter.headData = headData_;
	}

protected:
	FileStreamWriter* metaStream;
	FileStreamWriter* dataStream;

	void WriteFileHeader();
	void WriteFileFooter();
};


class ArchiveFileReader : public IArchiveFile
{
public:
	ArchiveFileReader();
	~ArchiveFileReader();

	void StartDecompress(const std::string& fileName_, ArchiveConfig& config_);
	bool ReadNextBin(DataChunk& buffer_, uint32& signature_);
	void FinishDecompress();


	// TODO: update accordingly to QVZ required data for decompression
	//
	const QualityCompressionData& GetQualityCompressionData() const
	{
		return fileFooter.quaData;
	}

	const FastqRawBlockStats::HeaderStats& GetHeadersCompressionData() const
	{
		return fileFooter.headData;
	}

protected:
	FileStreamReader* metaStream;
	FileStreamReader* dataStream;

	//uint64 blockIdx;
	std::vector<std::tuple<uint32, uint64, uint64> > blockArray;
	std::vector<std::tuple<uint32, uint64, uint64> >::const_iterator blockIterator;


	void ReadFileHeader();
	void ReadFileFooter();
};


#endif // H_DNARCHFILE
