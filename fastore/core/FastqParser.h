/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DNAPARSER
#define H_DNAPARSER

#include "Globals.h"

#include <vector>

#include "FastqRecord.h"
#include "Buffer.h"


/**
 * Parses records -- a general interface (deprecated!!!)
 *
 */
class IRecordParser
{
public:
	virtual ~IRecordParser() {}

	virtual bool ReadNextRecord(FastqRecord& rec_) = 0;
	virtual void WriteNextRecord(const FastqRecord& rec_) = 0;
};



/**
 * Parses DNA-only records (deprecated!!!)
 *
 */
class SingleDnaRecordParser : public IRecordParser
{
public:
	enum ParserMode
	{
		ParseRead,
		ParseWrite
	};

	SingleDnaRecordParser()
		:	memory(NULL)
		,	memoryPos(0)
		,	memorySize(0)
		,	skippedBytes(0)
		,	buf(NULL)
		,	stats(NULL)
	{}

	void StartParsing(DataChunk &chunk_, ParserMode mode_, FastqRawBlockStats* stats_ = NULL);

	virtual bool ReadNextRecord(FastqRecord& rec_);
	virtual void WriteNextRecord(const FastqRecord& rec_);

	uint64 FinishParsing(ParserMode mode_);
protected:
	byte* memory;
	uint64 memoryPos;
	uint64 memorySize;
	uint64 skippedBytes;

	Buffer* buf;
	FastqRawBlockStats* stats;

	bool ReadLine(uchar *str_, uint32& len_, uint32& size_);
	uint32 SkipLine();


	int32 Getc()
	{
		if (memoryPos == memorySize)
			return -1;
		return memory[memoryPos++];
	}

	void Skipc()
	{
		memoryPos++;
		skippedBytes++;
	}

	int32 Peekc()
	{
		if (memoryPos == memorySize)
			return -1;
		return memory[memoryPos];
	}

	void ExtendBuffer(uint64 minSize_)
	{
		ASSERT(buf != NULL);
		buf->Extend(MAX(minSize_, buf->Size() + (buf->Size() >> 1)), true);
		memory = buf->Pointer();
		memorySize = buf->Size();
	}
};


class SingleFastqRecordParser : public SingleDnaRecordParser
{
public:
	SingleFastqRecordParser(bool useHeaders_ = false, bool keepComments_ = false)
		:	useHeaders(useHeaders_)
		,	keepComments(keepComments_)
	{}

	bool ReadNextRecord(FastqRecord &rec_);
	void WriteNextRecord(const FastqRecord& rec_);

private:
	const bool useHeaders;
	const bool keepComments;
};



/**
 * Parses full FASTQ reads -- a general interface
 *
 */
class IRecordsParser
{
public:
	IRecordsParser(bool useHeaders_ = false, const std::string& libName_ = "SRX000000")
		:	useHeaders(useHeaders_)
		,	autoHeaderPrefix("@" + libName_ + ".")
	{}

	virtual ~IRecordsParser() {}

	virtual uint64 ParseTo(const std::vector<FastqRecord>& reads_, IFastqChunkCollection& chunk_, uint64 recStartIdx_ = 1) = 0;
	virtual uint64 ParseFrom(IFastqChunkCollection& chunks_, std::vector<FastqRecord>& records_, FastqRawBlockStats& stats_, bool keepComments_ = false) = 0;

	virtual void Reset() {}

protected:
	const bool useHeaders;
	const std::string autoHeaderPrefix;
};


class FastqRecordsParserSE : public IRecordsParser
{
public:
	using IRecordsParser::IRecordsParser;

	uint64 ParseTo(const std::vector<FastqRecord>& reads_, IFastqChunkCollection& chunk_, uint64 recStartIdx_ = 1);
	uint64 ParseFrom(IFastqChunkCollection& chunks_, std::vector<FastqRecord>& records_, FastqRawBlockStats& stats_, bool keepComments_ = false);
};


class FastqRecordsParserPE : public IRecordsParser
{
public:
	FastqRecordsParserPE(bool useHeaders_ = false, uint32 peFieldIdx_ = 0, const std::string& libName_ = "SRX000000")
		:	IRecordsParser(useHeaders_, libName_)
		,	peFieldIdx(peFieldIdx_)
	{}

	uint64 ParseTo(const std::vector<FastqRecord>& reads_, IFastqChunkCollection& chunk_, uint64 recStartIdx_ = 1);
	uint64 ParseFrom(IFastqChunkCollection& chunks_, std::vector<FastqRecord>& records_, FastqRawBlockStats& stats_, bool keepComments_ = false);

protected:
	const uint32 peFieldIdx;
};


/**
 * Dynmic parser stores the reads in a collection of multiple smaller chunks,
 * instead of a single one, large chunk of memory
 *
 */
class FastqRecordsParserDynSE : public FastqRecordsParserSE
{
public:
	using FastqRecordsParserSE::FastqRecordsParserSE;

	uint64 ParseTo(const std::vector<FastqRecord>& reads_, IFastqChunkCollection& chunk_, uint64 recStartIdx_ = 1);

private:
	static const uint64 DefaultBufferSize = 16 << 20;
};


class FastqRecordsParserDynPE : public FastqRecordsParserPE
{
public:
	using FastqRecordsParserPE::FastqRecordsParserPE;

	uint64 ParseTo(const std::vector<FastqRecord>& reads_, IFastqChunkCollection& chunk_, uint64 recStartIdx_ = 1);

private:
	static const uint64 DefaultBufferSize = 16 << 20;
};


#endif // H_DNAPARSER
