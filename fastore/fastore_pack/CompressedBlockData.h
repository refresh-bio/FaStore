/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_COMPRESSEDBLOCKDATA
#define H_COMPRESSEDBLOCKDATA

#include "../fastore_bin/Globals.h"
#include "../fastore_bin/Buffer.h"
#include "../fastore_bin/FastqRecord.h"


#include <vector>
#include <string>



/**
 * A compressed representation of binned reads
 * sharing the same signature
 *
 */
struct FastqCompressedBin
{
	static const uint64 DefaultBufferSize = 32 << 10;

	std::vector<DataChunk*> buffers;
	uint64 recordsCount;

	FastqCompressedBin(uint64 buffersNum_, uint64 bufferSize_ = DefaultBufferSize)
		:	recordsCount(0)
	{
		for (uint32 i = 0; i < buffersNum_; ++i)
			buffers.push_back(new DataChunk(bufferSize_));
	}

	virtual ~FastqCompressedBin()
	{
		for (uint32 i = 0; i < buffers.size(); ++i)
			delete buffers[i];
		buffers.clear();
	}

	virtual void Reset()
	{
		for (uint32 i = 0; i < buffers.size(); ++i)
			buffers[i]->size = 0;
		recordsCount = 0;

#if EXTRA_MEM_OPT
		for (uint32 i = 0; i < buffers.size(); ++i)
		{
			if (buffers[i]->data.Size() > DefaultBufferSize)
				buffers[i]->data.Shrink(DefaultBufferSize);
		}
#endif
	}
};



/**
 * A working buffers used for compressing
 *
 */
struct IFastqWorkBuffer
{
	FastqChunk fastqBuffer;
	FastqCompressedBin fastqWorkBin;

	IFastqWorkBuffer(uint32 buffersNum_, uint64 bufferSize_)
		:	fastqBuffer(bufferSize_)
		,	fastqWorkBin(buffersNum_)
	{}

	virtual ~IFastqWorkBuffer() {}

	void Reset()
	{
		fastqWorkBin.Reset();
		fastqBuffer.Reset();

#if EXTRA_MEM_OPT
		if (fastqBuffer.data.Size() > FastqChunk::DefaultBufferSize)
			fastqBuffer.data.Shrink(FastqChunk::DefaultBufferSize);
#endif
	}
};

struct FastqWorkBuffersSE : public IFastqWorkBuffer
{
	enum BufferNames
	{
		FlagBuffer = 0,
		LetterXBuffer,
		RevBuffer,

		HardReadsBuffer,
		LzIdBuffer,				// 5
		ShiftBuffer,
		MatchBuffer,
		MatchBinaryBuffer,

		TreeShiftBuffer,
		ConsensusMatchBuffer,	// 10
		ConsensusShiftBuffer,
		ConsensusLetterBuffer,

		QualityBuffer,			// ATM: quick solution

		ReadIdTokenBuffer,
		ReadIdValueBuffer,

		LastId = ReadIdValueBuffer
	};

	static const uint32 BuffersNum = LastId + 1;

	FastqWorkBuffersSE(uint32 buffersNum_ = BuffersNum,
					   uint64 bufferSize_ = FastqCompressedBin::DefaultBufferSize)
		:	IFastqWorkBuffer(buffersNum_, bufferSize_)
	{}

	static std::array<std::string, BuffersNum> GetBufferNames()
	{
		return {{"Flag", "LettersX", "Rev", "HardReads",
		 "LzId", "Shift", "Match", "MatchBinary", "TreeShift",
		 "CMatch", "CShift", "CLetters",
				"Quality", "ReadIdToken", "ReadIdValue"}};
	}
};



struct FastqWorkBuffersPE : public FastqWorkBuffersSE
{
	static const uint32 BuffersNum = FastqWorkBuffersSE::BuffersNum + 8;

	enum BufferNamesPE
	{
		FlagBufferPE = FastqWorkBuffersSE::LastId + 1,
		LetterXBufferPE,
		SwapFlagBuffer,
		HardReadsBufferPE,
		LzIdBufferPE,
		ShiftBufferPE,
		MatchRLEBufferPE,
		MatchBinaryBufferPE
	};

	FastqWorkBuffersPE(uint32 buffersNum_ = BuffersNum,
					   uint64 bufferSize_ = FastqCompressedBin::DefaultBufferSize)
		:	FastqWorkBuffersSE(buffersNum_, bufferSize_)
	{}

	static std::array<std::string, BuffersNum> GetBufferNames()
	{
		return {{"Flag", "LettersX", "Rev", "HardReads",
		 "LzId", "Shift", "Match", "MatchBinary", "TreeShift",
		 "CMatch", "CShift", "CLetters", "Quality", "ReadIdToken", "ReadIdValue",
		 "PE_Flag", "PE_LettersX", "PE_Swap", "PE_Hard",
		 "PE_LzId", "PE_Shift", "PE_MatchRLE", "PE_MatchBinary"}};
	}
};


struct CompressedFastqBlockStats
{
	static const uint32 BuffersNum = FastqWorkBuffersSE::BuffersNum;

	// TODO: optimize using enums
	//
	std::map<std::string, uint64> counts;
	std::map<std::string, std::map<int32, uint64> > freqs;
	std::map<std::string, std::vector<uint64> > bufferSizes;

	uint64 recordsCount;
	uint32 currentSignature;

	CompressedFastqBlockStats()
		:	recordsCount(0)
		,	currentSignature(0)
    {}

    void Reset()
    {
		counts.clear();
		freqs.clear();
		bufferSizes.clear();

		recordsCount = 0;
		currentSignature = 0;
    }

	void Update(const CompressedFastqBlockStats& stats_)
	{
		for (auto sv : stats_.counts)
			counts[sv.first] += sv.second;

		for (const auto& sv : stats_.freqs)
		{
			auto& ff = freqs[sv.first];
			for (auto nf : sv.second)
				ff[nf.first] += nf.second;
		}

		for (const auto& sv : stats_.bufferSizes)
		{
			auto& buf = bufferSizes[sv.first];
			if (buf.size() != sv.second.size())
				buf.resize(sv.second.size(), 0);

			for (uint32 i = 0; i < buf.size(); ++i)
				buf[i] += sv.second[i];
		}
	}
};




/**
 * A compressed block of FASTQ reads keeping only raw binary data
 *
 */
struct CompressedFastqBlock
{
	uint32 signatureId;
	FastqChunk dataBuffer;

	std::string log;

	CompressedFastqBlockStats stats;

	CompressedFastqBlock(uint64 bufferSize_ = FastqChunk::DefaultBufferSize)
		:	signatureId(0)
		,	dataBuffer(bufferSize_)
	{}

	void Reset()
	{
		signatureId = 0;
        stats.Reset();
		dataBuffer.Reset();

		log.clear();

#if EXTRA_MEM_OPT
		if (dataBuffer.data.Size() > FastqChunk::DefaultBufferSize)
			dataBuffer.data.Shrink(FastqChunk::DefaultBufferSize);
#endif
	}
};


#endif // H_COMPRESSEDBLOCKDATA
