/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef FASTQRECORD_H
#define FASTQRECORD_H

#include "Globals.h"
#include "Buffer.h"

#include <vector>
#include <map>
#include <string.h>
#include <array>


/**
 * A definition of FASTQ record with sample operations
 * WARN: it does not store the DNA, ID or QUA data, it only
 * links to the memory containig it
 *
 */
struct FastqRecord
{
	static const uint32 MaxSeqLen = 256;			// temporary, can be increased
	static const uint32 MaxTagLen = 255;
	static const uint32 DefaultQualityOffset = 33;

	enum RecordFlags
	{
		FlagReadIsReverse	= BIT(0),
		FlagIsPairSwapped	= BIT(1),
	};


	char* seq;
	char* qua;
	char* head;

	uint16 seqLen;
	uint16 auxLen;		// for PE
	uint16 minimPos;
	uint8 headLen;
	uint8 flags;

	FastqRecord()
		:	seq(NULL)
		,	qua(NULL)
		,	head(NULL)
		,	seqLen(0)
		,	auxLen(0)
		,	minimPos(0)
		,	headLen(0)
		,	flags(0)
	{}

	static const char* GetRCCodes()
	{

#ifdef __APPLE__
		static const char
#else
		static const thread_local char
#endif
			rcCodes[24] = {	-1,'T',-1,'G',-1,-1,-1,'C',		// 64+
							-1,-1,-1,-1,-1,-1,'N',-1,		// 72+,
							-1,-1,-1,-1,'A',-1,-1,-1,		// 80+
						};

		return rcCodes;
	}

	// TODO: deprecated
	//
	void ComputeRC(FastqRecord& rc_) const
	{
		const char* rcCodes = FastqRecord::GetRCCodes();

		ASSERT(rc_.seq != NULL);

		ASSERT(seqLen + auxLen > 0);
		ASSERT(auxLen == 0 || seqLen == auxLen);		// so far only same length reads are accepted

		rc_.seqLen = (auxLen > 0) ? auxLen : seqLen;
		rc_.auxLen = (auxLen > 0) ? seqLen : auxLen;

		const uint32 len = seqLen + auxLen;
		for (uint32 i = 0; i < len; ++i)				// auxLen to handle paired-end
		{
			ASSERT(seq[i] > 64 && seq[i] < 88);
			rc_.seq[len-1-i] = rcCodes[(int32)seq[i] - 64];
			ASSERT(rc_.seq[len-1-i] != -1);
		}

		if (qua != NULL)
		{
			for (uint32 i = 0; i < len; ++i)
				rc_.qua[len-1-i] = qua[i];
		}
		else
		{
			rc_.qua = NULL;
		}

		rc_.minimPos = minimPos;
	}

	void Reset()
	{
		flags = 0;
		minimPos = 0;
	}

	// TODO:
	// void Revert();

        void CopyFrom(const FastqRecord& rec_, bool copyHeader_ = false)
	{
		ASSERT(rec_.seqLen == seqLen);
		ASSERT(rec_.auxLen == auxLen);
		std::copy(rec_.seq, rec_.seq + rec_.seqLen + rec_.auxLen, seq);

		if (rec_.qua != NULL)
		{
			ASSERT(qua != NULL);
			std::copy(rec_.qua, rec_.qua + rec_.seqLen + rec_.auxLen, qua);
		}

                if (rec_.head != NULL && copyHeader_)
                {
                    ASSERT(head != NULL);
                    std::copy(rec_.head, rec_.head + rec_.headLen, head);
                    headLen = rec_.headLen;
                }
	}


	// flags handling
	//
	bool IsSetFlag(uint32 flag_) const
	{
		return (flags & flag_) != 0;
	}

	void SetFlag(uint32 flag_, bool b_)
	{
		b_ ? (flags |= flag_) : (flags &= ~flag_);
	}

	bool IsReadReverse() const
	{
		return IsSetFlag(FlagReadIsReverse);
	}

	void SetReadReverse(bool state_)
	{
		SetFlag(FlagReadIsReverse, state_);
	}

	bool IsPairSwapped() const
	{
		return IsSetFlag(FlagIsPairSwapped);
	}

	void SetPairSwapped(bool state_)
	{
		SetFlag(FlagIsPairSwapped, state_);
	}


	// PE-related methods
	//
	FastqRecord GetPair() const
	{
		FastqRecord rec;
		rec.seq = seq + seqLen;
		rec.qua = qua + seqLen;
		rec.seqLen = auxLen;
		rec.head = head;
		rec.headLen = headLen;

		return rec;
	}

	void SwapReads()
	{
		ASSERT(seqLen == auxLen);

		std::swap_ranges(seq, seq + seqLen, seq + seqLen);
		std::swap_ranges(qua, qua + seqLen, qua + seqLen);


		SetPairSwapped(!IsPairSwapped());
	}
};


/**
 * Buffers the SEQ and DNA information of a FASTQ read
 *
 */
struct FastqRecordBuffer : public FastqRecord
{
	FastqRecordBuffer()
	{
		seq = seqBuffer.data();
		qua = quaBuffer.data();
	}

	std::array<char, MaxSeqLen * 2> seqBuffer;
	std::array<char, MaxSeqLen * 2> quaBuffer;
};


/**
 * Comparator used for sorting FASTQ reads, when performin matching
 *
 */
struct IFastqComparator
{
	bool CompareReads(const FastqRecord& r1_, const FastqRecord& r2_)
	{
		const char* r1p = r1_.seq + r1_.minimPos;
		const char* r2p = r2_.seq + r2_.minimPos;

		// WARN: quick'n'dirty PE fix
		uint32 len = MIN(r1_.seqLen + r1_.auxLen - r1_.minimPos, r2_.seqLen + r2_.auxLen - r2_.minimPos);

		const int32 r = strncmp(r1p, r2p, len);
		if (r == 0)
		{
			// we can extend this one
			ASSERT(r1_.minimPos < r1_.seqLen);
			ASSERT(r2_.minimPos < r2_.seqLen);

			// reverse-compare for better catching exact matches
			if (r1_.minimPos == r2_.minimPos)
			{
				for (int32 i = (int32)r1_.minimPos; i >= 0; i--)
				{
					if (r1_.seq[i] < r2_.seq[i])
						return true;
					if (r1_.seq[i] > r2_.seq[i])
						return false;
				}
				return false;	// equal
			}

			return r1_.minimPos > r2_.minimPos;
		}
		return r < 0;
	}
};

template <typename _TRecordType>
struct TFastqComparator
{};

template <>
struct TFastqComparator<const FastqRecord*> : public IFastqComparator
{
	bool operator() (const FastqRecord* r1_, const FastqRecord* r2_)
	{
		return CompareReads(*r1_, *r2_);
	}
};

template <>
struct TFastqComparator<const FastqRecord&> : public IFastqComparator
{
	bool operator() (const FastqRecord& r1_, const FastqRecord& r2_)
	{
		return CompareReads(r1_, r2_);
	}
};

typedef TFastqComparator<const FastqRecord&> FastqComparator;
typedef TFastqComparator<const FastqRecord*> FastqComparatorPtr;



/**
 * Stores the basic statistics when binning FASTQ records
 */
struct FastqRecordBinStats
{
	uint32 minSeqLen;
	uint32 maxSeqLen;

	uint32 minAuxLen;
	uint32 maxAuxLen;

	// here we can also add extra dna / quality stats which will be calculated per-bin
	//

	FastqRecordBinStats()
	{
		Clear();
	}

	void Clear()
	{
		minSeqLen = minAuxLen = (uint32)-1;
		maxSeqLen = maxAuxLen = 0;
	}

	void Update(const FastqRecord& rec_)
	{
		maxSeqLen = MAX(maxSeqLen, rec_.seqLen);
		minSeqLen = MIN(minSeqLen, rec_.seqLen);

		maxAuxLen = MAX(maxAuxLen, rec_.auxLen);
		minAuxLen = MIN(minAuxLen, rec_.auxLen);
	}

	void Update(const FastqRecordBinStats& stats_)
	{
		maxSeqLen = MAX(maxSeqLen, stats_.maxSeqLen);
		minSeqLen = MIN(minSeqLen, stats_.maxSeqLen);

		maxAuxLen = MAX(maxAuxLen, stats_.maxAuxLen);
		minAuxLen = MIN(minAuxLen, stats_.minAuxLen);
	}
};


/**
 * A collection of memory chunks for storing FASTQ data
 *
 */
typedef DataChunk FastqChunk;

struct IFastqChunkCollection
{
	IFastqChunkCollection(uint32 chunkNum_ = 0, uint64 bufferSize_ = 0)
		:	defaultBufferSize(bufferSize_)
	{
		ASSERT(bufferSize_ != 0 || chunkNum_ == 0);
		for (uint32 i = 0; i < chunkNum_; ++i)
			chunks.push_back(new FastqChunk(defaultBufferSize));
	}

	virtual ~IFastqChunkCollection()
	{
		Clear();
	}

	void Reset()
	{
		for (FastqChunk* c : chunks)
			c->Reset();
	}

	void Clear()
	{
		for (FastqChunk* c : chunks)
			delete c;
		chunks.clear();
	}

	const uint64 defaultBufferSize;
	std::vector<FastqChunk*> chunks;
};


struct FastqChunkCollectionSE : public IFastqChunkCollection
{
	FastqChunkCollectionSE(uint64 bufferSize_ = FastqChunk::DefaultBufferSize)
		:	IFastqChunkCollection(1, bufferSize_)
	{}
};


struct FastqChunkCollectionPE : public IFastqChunkCollection
{
	enum ChunkNames
	{
		InputChunk1 = 0,
		InputChunk2,
		OutputChunk
	};

	FastqChunkCollectionPE(uint64 bufferSize_ = FastqChunk::DefaultBufferSize)
		:	IFastqChunkCollection(3, bufferSize_)
	{}
};


/**
 * Represents a single bin containing pointers to the FASTQ records
 * falling into it
 *
 */
struct FastqRecordsPtrBin
{
	std::vector<FastqRecord*> records;
	FastqRecordBinStats stats;

	void Clear()
	{
		records.clear();
		stats.Clear();

#if EXTRA_MEM_OPT
		records.shrink_to_fit();
#endif
	}
};


/**
 * Data types used to represent the relation graph between
 * the matched FASTQ reads
 *
 */
struct MatchNode;

struct MatchNodesPtrBin
{
	std::vector<const MatchNode*> nodes;
	FastqRecordBinStats stats;				// do we need this one when the lenght is fixed?

	void Clear()
	{
		nodes.clear();
		stats.Clear();

#if EXTRA_MEM_OPT
		nodes.shrink_to_fit();
#endif
	}
};


#endif // FASTQRECORD_H
