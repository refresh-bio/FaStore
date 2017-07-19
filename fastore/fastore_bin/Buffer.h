/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BUFFER
#define H_BUFFER

#include "Globals.h"

#include <algorithm>

#include "Utils.h"


/**
 * Memory buffer interface
 *
 */
class IBuffer
{
public:
	IBuffer(byte* ptr_, uint64 size_)
		:	buffer(ptr_)
		,	size(size_)
	{}

	virtual ~IBuffer() {}

	uint64 Size() const
	{
		return size;
	}

	byte* Pointer() const
	{
		return (byte*)buffer;
	}

	uint64* Pointer64() const
	{
		return (uint64*)buffer;
	}

	void SetSize(uint64 s_)
	{
		size = s_;
	}

	void SetPointer(byte* ptr_)
	{
		buffer = ptr_;
	}

protected:
	byte* buffer;
	uint64 size;

private:
	IBuffer(const IBuffer& )
	{}

	IBuffer& operator= (const IBuffer& )
	{ return *this; }
};


/**
 * Memory buffer implementing RAII-style memory management
 *
 */
class Buffer : public IBuffer
{
public:
	Buffer(uint64 size_)
		:	IBuffer(Alloc(size_), size_)
	{
		ASSERT(size_ != 0);
	}

	~Buffer()
	{
		delete[] buffer;
	}

	void Extend(uint64 size_, bool copy_ = false)
	{
		ASSERT(size < size_);

                byte* p = new (std::nothrow) byte[size_];
                ASSERT(p != NULL);

		if (copy_)
			std::copy(buffer, buffer + size, (byte*)p);

		delete[] buffer;

		buffer = (byte*)p;
		size = size_;
	}

	void Shrink(uint64 size_)
	{
		byte* p = new byte[size_];
		delete[] buffer;

		buffer = (byte*)p;
		size = size_;
	}

	void Swap(Buffer& b)
	{
		std::swap(b.buffer, buffer);
		std::swap(b.size, size);
	}

	static byte* Alloc(uint64 size_)
	{
		uint64 size64 = size_ / 8;
		if (size64 * 8 < size_)
			size64 += 1;
		return (byte*)(new uint64[size64]);
	}

private:
	Buffer(const Buffer& ) : IBuffer(NULL, 0)
	{}

	Buffer& operator= (const Buffer& )
	{ return *this; }
};



/**
 * Keeps memory buffer with the size of the used data
 *
 */
struct DataChunk
{
	static const uint64 DefaultBufferSize = 256 << 10;

	Buffer data;
	uint64 size;

	DataChunk(uint64 bufferSize_ = DefaultBufferSize)
		:	data(bufferSize_)
		,	size(0)
	{}

	void Reset()
	{
		size = 0;
	}
};


/**
 * RAII-style array of data chunks
 *
 */
struct DataChunkCollection
{
	std::vector<DataChunk*> chunks;

	DataChunkCollection()
	{}

	~DataChunkCollection()
	{
		for (auto dc : chunks)
			delete dc;
	}

	DataChunk* AddNewChunk(uint64 chunkSize_ = DataChunk::DefaultBufferSize)
	{
		DataChunk* dc = new DataChunk(chunkSize_);
		ASSERT(dc != NULL);
		chunks.push_back(dc);
		return dc;
	}
};


#endif // H_BUFFER
