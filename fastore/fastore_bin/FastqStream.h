/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_FASTQSTREAM
#define H_FASTQSTREAM

#include "Globals.h"

#include <string>
#include <vector>

#include "FileStream.h"
#include "Exception.h"
#include "FastqRecord.h"


/**
 * Reads FASTQ file(s) chunk-wise -- a general interface
 *
 */
class IFastqStreamReaderBase
{
public:
	IFastqStreamReaderBase()
		:	usesCrlf(false)
	{}

	virtual ~IFastqStreamReaderBase()
	{}

	virtual bool ReadNextChunk(IFastqChunkCollection& chunk_) = 0;

	virtual bool Eof() const = 0;
	virtual void Close() = 0;

protected:
	bool usesCrlf;

	uint64 GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_);

	void SkipToEol(uchar* data_, uint64& pos_, const uint64 size_)
	{
		ASSERT(pos_ < size_);

		while (data_[pos_] != '\n' && data_[pos_] != '\r' && pos_ < size_)
			++pos_;

		if (data_[pos_] == '\r' && pos_ < size_)
		{
			if (data_[pos_ + 1] == '\n')
			{
				usesCrlf = true;
				++pos_;
			}
		}
	}

	int64 Read(IDataStreamReader* stream_, byte* memory_, uint64 size_)
	{
		ASSERT(stream_ != NULL);
		ASSERT(memory_ != NULL);
		return stream_->Read(memory_, size_);
	}
};


class IFastqStreamReaderSE : public IFastqStreamReaderBase
{
public:
	IFastqStreamReaderSE(uint64 maxReadBufferSize_ = MaxReadBufferSize)
		:	maxReadBufferSize(maxReadBufferSize_)
		,	stream(NULL)
		,	readBuffer(maxReadBufferSize_)
		,	readBufferSize(0)
		,	eof(false)
	{}

	bool Eof() const
	{
		return eof;
	}

	bool ReadNextChunk(IFastqChunkCollection& chunk_);

	void Close()
	{
		ASSERT(stream != NULL);
		stream->Close();
	}


protected:
	static const uint32 MaxReadBufferSize = 1 << 13;

	const uint64 maxReadBufferSize;

	IDataStreamReader* stream;
	Buffer readBuffer;
	uint64 readBufferSize;
	bool eof;

	int64 Read(byte* memory_, uint64 size_)
	{
		return IFastqStreamReaderBase::Read(stream, memory_, size_);
	}
};


class IFastqStreamReaderPE : public IFastqStreamReaderSE
{
public:
	IFastqStreamReaderPE(uint64 maxReadBufferSize_ = MaxPairBufferSize)
		:	IFastqStreamReaderSE(maxReadBufferSize_)
		,	stream_2(NULL)
		,	pairBuffer(maxReadBufferSize_)
		,	pairBufferSize(0)
		,	eof_2(false)
	{}

	bool Eof() const
	{
		return IFastqStreamReaderSE::Eof() && eof_2;
	}

	bool ReadNextChunk(IFastqChunkCollection& chunk_);

	void Close()
	{
		IFastqStreamReaderSE::Close();
		ASSERT(stream_2 != NULL);
		stream_2->Close();
	}


protected:
	static const uint32 MaxPairBufferSize = 1 << 20;		// TODO: take as fraction of a MAX input buffer size

	IDataStreamReader* stream_2;
	Buffer pairBuffer;
	uint64 pairBufferSize;
	bool eof_2;

	int64 Read_1(byte* memory_, uint64 size_)
	{
		return IFastqStreamReaderSE::Read(memory_, size_);
	}

	int64 Read_2(byte* memory_, uint64 size_)
	{
		return IFastqStreamReaderBase::Read(stream_2, memory_, size_);
	}

	uint64 ParseNextReadId(const uchar* data_, uint64 maxLen_);

private:
	using IFastqStreamReaderSE::ReadNextChunk;
};



/**
 * Writes FASTQ file(s) chunk-wise -- a general interface
 *
 */
class IFastqStreamWriter
{
public:
	virtual ~IFastqStreamWriter()
	{}

	virtual void WriteNextChunk(const IFastqChunkCollection& chunk_) = 0;

	virtual void Close() = 0;

protected:
	int64 Write(IDataStreamWriter* stream_, const DataChunk* chunk_)
	{
		ASSERT(chunk_ != NULL);
		ASSERT(stream_ != NULL);
		return stream_->Write(chunk_->data.Pointer(), chunk_->size);
	}
};


class IFastqStreamWriterSE : public IFastqStreamWriter
{
public:
	IFastqStreamWriterSE()
		:	stream(NULL)
	{}

	void WriteNextChunk(const IFastqChunkCollection& chunk_)
	{
		ASSERT(chunk_.chunks.size() >= 1);

		for (DataChunk* c : chunk_.chunks)
		{
			// highly probable that when we reach first zero-size chunk,
			// the restil will be empty too
			if (c->size > 0)
				Write(stream, c);
		}
	}

	void Close()
	{
		ASSERT(stream != NULL);
		stream->Close();
	}

protected:
	IDataStreamWriter* stream;
};


class IFastqStreamWriterPE : public IFastqStreamWriterSE
{
public:
	IFastqStreamWriterPE()
		:	stream_2(NULL)
	{}

	void WriteNextChunk(const IFastqChunkCollection& chunk_)
	{
		ASSERT(chunk_.chunks.size() >= 2);

		if (chunk_.chunks.size() == 3) // special case for binning / rebinning
		{
			ASSERT(chunk_.chunks[0]->size > 0);
			ASSERT(chunk_.chunks[1]->size > 0);
			Write(stream, chunk_.chunks[0]);
			Write(stream_2, chunk_.chunks[1]);

			return;
		}

		uint32 i = 0;
		while (i < chunk_.chunks.size())
		{
			DataChunk* dc1 = chunk_.chunks[i];
			if (dc1->size > 0)
				Write(stream, dc1);

			ASSERT(chunk_.chunks.size() > i + 1);
			DataChunk* dc2 = chunk_.chunks[i+1];
			if (dc2->size > 0)
				Write(stream_2, dc2);

			i += 2;
		}
	}

	void Close()
	{
		IFastqStreamWriterSE::Close();

		ASSERT(stream_2 != NULL);
		stream_2->Close();
	}


protected:
	IDataStreamWriter* stream_2;


private:
	// hide:
	using IFastqStreamWriterSE::WriteNextChunk;
};



/**
 * Wrappers over FASTQ reader(s)/writers(s)
 *
 */
template <class _TStreamInterface, class _TStream, class _TStreamInput>
class TFastqStreamSE : public _TStreamInterface
{
public:
	TFastqStreamSE(const _TStreamInput& input_)
	{
		_TStreamInterface::stream = new _TStream(input_);
	}

	~TFastqStreamSE()
	{
		delete _TStreamInterface::stream;
	}
};


template <class _TStreamInterface, class _TStream, class _TStreamInput>
class TFastqStreamPE : public _TStreamInterface
{
public:
	TFastqStreamPE(const _TStreamInput& input1_, const _TStreamInput& input2_)
	{

		_TStreamInterface::stream = new _TStream(input1_);
		try
		{
			_TStreamInterface::stream_2 = new _TStream(input2_);
		}
		catch (const Exception& e_)
		{
			delete _TStreamInterface::stream;
			_TStreamInterface::stream = NULL;
			throw e_;
		}
	}

	~TFastqStreamPE()
	{
		delete _TStreamInterface::stream;
		delete _TStreamInterface::stream_2;
	}
};


// single FASTQ file reading/writing both SE and PE
//
typedef TFastqStreamSE<IFastqStreamReaderSE, FileStreamReader, std::string> FastqFileReaderSE;
typedef TFastqStreamSE<IFastqStreamWriterSE, FileStreamWriter, std::string> FastqFileWriterSE;

typedef TFastqStreamPE<IFastqStreamReaderPE, FileStreamReader, std::string> FastqFileReaderPE;
typedef TFastqStreamPE<IFastqStreamWriterPE, FileStreamWriter, std::string> FastqFileWriterPE;


// multi FASTQ file reading both SE and PE for both raw and gz-compressed
//
typedef TFastqStreamSE<IFastqStreamReaderSE, MultiFileStreamReader, std::vector<std::string> > MultiFastqFileReaderSE;
typedef TFastqStreamPE<IFastqStreamReaderPE, MultiFileStreamReader, std::vector<std::string> > MultiFastqFileReaderPE;

typedef TFastqStreamSE<IFastqStreamReaderSE, MultiFileStreamReaderGz, std::vector<std::string> > MultiFastqFileReaderGzSE;
typedef TFastqStreamPE<IFastqStreamReaderPE, MultiFileStreamReaderGz, std::vector<std::string> > MultiFastqFileReaderGzPE;


#endif // H_FASTQSTREAM
