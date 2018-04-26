/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINFILE
#define H_BINFILE

#include "../core/Globals.h"

#include <vector>

#include "Params.h"

#include "../core/BinBlockData.h"
#include "../core/FileStream.h"
#include "../qvz/QVZ.h"


/**
 * Interface for reading/writing files containing bined data
 *
 */
class IBinFile
{
public:
	struct BinFileFooter
	{
		static const uint64 ParametersSize = sizeof(BinModuleConfig);

		struct BlockMetaData : public BinaryBinDescriptor
		{
			uint64 metaFileOffset;
			uint64 dnaFileOffset;
			uint64 quaFileOffset;
			uint64 headFileOffset;							// optional

			BlockMetaData()
				:	metaFileOffset(0)
				,	dnaFileOffset(0)
				,	quaFileOffset(0)
				,	headFileOffset(0)
			{}

			BlockMetaData(const BinaryBinDescriptor& b_)
				:	BinaryBinDescriptor(b_)
				,	metaFileOffset(0)
				,	dnaFileOffset(0)
				,	quaFileOffset(0)
				,	headFileOffset(0)
			{}
		};

		struct BinInfo
		{
			std::vector<BlockMetaData> blocksMetaData;

			uint64 totalMetaSize;
			uint64 totalDnaSize;
			uint64 totalQuaSize;
			uint64 totalHeadSize;							// optional

			uint64 totalRawDnaSize;
			uint64 totalRawHeadSize;						// optional
			uint64 totalRecordsCount;

			BinInfo()
				:	totalMetaSize(0)
				,	totalDnaSize(0)
				,	totalQuaSize(0)
				,	totalHeadSize(0)
				,	totalRawDnaSize(0)
				,	totalRawHeadSize(0)
				,	totalRecordsCount(0)
			{}
		};

		BinModuleConfig params;

		std::map<uint32, BinInfo> binOffsets;

		QualityCompressionData quaData;

		// TODO: header compression data
		//
		FastqRawBlockStats::HeaderStats headData;


		void Clear()
		{
			binOffsets.clear();
		}
	};

	IBinFile()
	{
		std::fill((uchar*)&fileHeader, (uchar*)&fileHeader + sizeof(BinFileHeader), 0);
	}

	virtual ~IBinFile() {}

protected:
	struct BinFileHeader
	{
		static const uint64 ReservedBytes = 7;
		static const uint64 HeaderSize = 4*8 + 1 + ReservedBytes;

		uint64 footerOffset;
		uint64 recordsCount;
		uint64 blockCount;
		uint64 footerSize;
		bool usesHeaderStream;

		uchar reserved[ReservedBytes];
	};

	BinFileHeader fileHeader;
	BinFileFooter fileFooter;
};



class BinFileWriter : public IBinFile
{
public:
	BinFileWriter();
	~BinFileWriter();

	void StartCompress(const std::string& filename_,
					   const BinModuleConfig& params_);

	void WriteNextBlock(const BinaryBinBlock* block_);
	void FinishCompress();

	const BinFileFooter& GetFileFooter() const
	{
		return fileFooter;
	}

	const BinFileHeader& GetFileHeader() const
	{
		return fileHeader;
	}

	void SetQualityCompressionData(const QualityCompressionData& qua_)
	{
		fileFooter.quaData = qua_;
	}

	void SetHeaderCompressionData(const FastqRawBlockStats::HeaderStats& head_)
	{
		fileFooter.headData = head_;
	}

protected:
	IDataStreamWriter* metaStream;
	IDataStreamWriter* dnaStream;
	IDataStreamWriter* quaStream;
	IDataStreamWriter* headStream;

	FastqRawBlockStats globalFastqStats;

	void WriteFileHeader();
	void WriteFileFooter();
};



class BinFileReader : public IBinFile
{
public:
	BinFileReader();
	~BinFileReader();

	void StartDecompress(const std::string& fileName_, BinModuleConfig& params_);

	bool ReadNextBlock(BinaryBinBlock* block_);
	void FinishDecompress();

	uint64 BlockCount() const
	{
		return fileHeader.blockCount;
	}

	const BinFileFooter& GetFileFooter() const
	{
		return fileFooter;
	}

protected:
	IDataStreamReader* metaStream;
	IDataStreamReader* dnaStream;
	IDataStreamReader* quaStream;
	IDataStreamReader* headStream;

	std::map<uint32, BinFileFooter::BinInfo>::const_iterator metaOffsetIterator;

	void ReadBlock(uint32 signature_, BinaryBinBlock* block_);

	void ReadFileHeader();
	void ReadFileFooter();
};

#endif // H_BINFILE
