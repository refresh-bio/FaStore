/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BITMEMORY
#define H_BITMEMORY

#include "Globals.h"

#include <vector>
#include <algorithm>

#include "Buffer.h"
#include "Utils.h"



/**
 * Reads the data from a specified memory buffer
 *
 */
class BitMemoryReader
{
public:
	BitMemoryReader(const Buffer& buffer_, uint64 size_, uint64 offset_ = 0)
		:	buffer(buffer_)
		,	size(size_)
		,	position(offset_)
		,	wordBuffer(0)
		,	wordBufferPos(0)
	{}

	uint64 Size() const
	{
		return size;
	}

	uint64 Position() const
	{
		return position;
	}

	void SetPosition(uint64 pos_)
	{
		ASSERT(pos_ <= size);
		position = pos_;
	}

	byte* Pointer() const
	{
		return buffer.Pointer();
	}

	// TODO: refactor
	uint32 GetBit()
	{
		if (wordBufferPos == 0)
		{
			wordBuffer = GetByte();
			wordBufferPos = 7;
			return (wordBuffer >> 7) & 1;
		}

		return (wordBuffer >> (--wordBufferPos)) & 1;
	}

	uint32 Get2Bits()
	{
		if (wordBufferPos >= 2)
		{
			wordBufferPos -= 2;
			return (wordBuffer >> wordBufferPos) & 3;
		}

		if (wordBufferPos == 0)
		{
			wordBuffer = GetByte();
			wordBufferPos = 6;
			return (wordBuffer >> wordBufferPos) & 3;
		}

		uint32 word = (wordBuffer & 1) << 1;
		wordBuffer = GetByte();
		wordBufferPos = 7;
		word += wordBuffer >> wordBufferPos;
		return word & 3;
	}

	uint32 GetBits(uint32 n_)
	{
		ASSERT(n_ > 0 && n_ < 32);

		uint32 word = 0;
		while (n_)
		{
			if (wordBufferPos == 0)
			{
				wordBuffer = GetByte();
				wordBufferPos = 8;
			}

			if (n_ > wordBufferPos)
			{
				word <<= wordBufferPos;
				word += wordBuffer & BitMask(wordBufferPos);
				n_ -= wordBufferPos;
				wordBufferPos = 0;
			}
			else
			{
				word <<= n_;
				wordBufferPos -= n_;
				word += (wordBuffer >> wordBufferPos) & BitMask(n_);
				break;
			}
		}
		return word;
	}

	byte GetByte()
	{
		ASSERT(position < size);
		return (buffer.Pointer())[position++];
	}

	int32 PeekByte() const
	{
		if (position < size)
			return (buffer.Pointer())[position];
		return -1;
	}

	void GetBytes(uchar *data, uint64 n_bytes)
	{
		ASSERT(position + n_bytes <= size);

		std::copy(buffer.Pointer() + position, buffer.Pointer() + position + n_bytes, data);
		position += n_bytes;
	}

	uint32 Get4Bytes()
	{
		uint32 c = GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		return (c << 8) | GetByte();
	}

	uint64 Get8Bytes()
	{
		uint64 c = GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		c = (c << 8) | GetByte();
		return c;
	}

	uint32 Get2Bytes()
	{
		return GetByte() << 8 | GetByte();
	}

	void SkipBytes(uint64 n_)
	{
		ASSERT(position + n_ < size);

		position += n_;
	}

	void FlushInputWordBuffer()
	{
		wordBufferPos = 0;
	}

	void Reset(uint64 size_ = 0)
	{
		ASSERT(buffer.Size() >= size_);

		size = size_;
		position = 0;
		wordBufferPos = 0;
		wordBuffer = 0;
	}

private:
	static const uint32 WordBufferSize = 8;

	const Buffer& buffer;
	uint64 size;
	uint64 position;

	uint32 wordBuffer;
	uint32 wordBufferPos;

	uint32 BitMask(uint32 n_)
	{
		ASSERT(n_ < 32);
		return ((uint32)1 << n_) - 1;
	}
};


/**
 * Writes the data to a specified memory buffer
 *
 */
class BitMemoryWriter
{
public:
	static const uint32 DefaultBufferSize = 1 << 20;

	BitMemoryWriter(Buffer& buffer_, uint64 offset_ = 0)
		:	buffer(buffer_)
		,	position(offset_)
		,	wordBuffer(0)
		,	wordBufferPos(0)
	{
		ASSERT(position <= buffer.Size());
	}

	~BitMemoryWriter()
	{
	}

	uint64 Position() const
	{
		return position;
	}

	void SetPosition(uint64 pos_)
	{
		ASSERT(pos_ <= buffer.Size());
		position = pos_;
	}

	byte* Pointer() const
	{
		return buffer.Pointer();
	}

	template <typename _T>
	void PutBit(_T b_)
	{
		if (wordBufferPos < WordBufferSize)
		{
			wordBuffer <<= 1;
			wordBuffer += ((uint32)b_) & 1;
			++wordBufferPos;
		}
		else
		{
			Put4Bytes(wordBuffer);
			wordBufferPos = 1;
			wordBuffer = ((uint32)b_) & 1;
		}
	}

	void Put2Bits(uint32 word_)
	{
		ASSERT(word_ <= 3);
		word_ &= 3;

		if (wordBufferPos + 2 <= WordBufferSize)
		{
			wordBuffer <<= 2;
			wordBuffer += word_;
			wordBufferPos += 2;
		}
		else if (wordBufferPos == WordBufferSize)
		{
			Put4Bytes(wordBuffer);
			wordBufferPos = 2;
			wordBuffer = word_;
		}
		else
		{
			wordBuffer <<= 1;
			wordBuffer += word_ >> 1;
			Put4Bytes(wordBuffer);
			wordBuffer = word_ & 1;
			wordBufferPos = 1;
		}
	}

	void PutBits(uint32 word_, uint32 n_)
	{
		ASSERT(n_ > 0);
		word_ &= BitMask(n_);

		uint32 rest = WordBufferSize - wordBufferPos;
		if (n_ >= rest)
		{
			n_ -= rest;
			wordBuffer <<= rest;
			wordBuffer += word_ >> n_;
			wordBufferPos = 0;
			Put4Bytes(wordBuffer);
			wordBuffer = 0;
		}

		wordBuffer <<= n_;
		wordBuffer += word_ & BitMask(n_);
		wordBufferPos += n_;
	}

	void PutByte(byte b_)
	{
		if (position >= buffer.Size())
			buffer.Extend(buffer.Size() + (buffer.Size() >> 2), true);

		(buffer.Pointer())[position++] = b_;
	}

	void Put2Bytes(uint32 word_)
	{
		PutByte((uchar)(word_ >> 8));
		PutByte((uchar)(word_ & 0xFF));
	}

	void PutBytes(const byte *data_, uint64 n_)
	{
		if (position + n_ > buffer.Size())
			buffer.Extend(position + n_ + (buffer.Size() >> 2), true);

		std::copy(data_, data_ + n_, buffer.Pointer() + position);
		position += n_;
	}

	void FillBytes(byte data_, uint64 n_)
	{
		if (position + n_ > buffer.Size())
			buffer.Extend(position + n_ + (buffer.Size() >> 2), true);

		std::fill(buffer.Pointer() + position, buffer.Pointer() + position + n_, data_);
		position += n_;
	}

	void Put8Bytes(uint64 data_)
	{
		PutByte(data_ >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
		PutByte((data_ <<= 8) >> 56);
	}

	void Put4Bytes(uint32 data_)
	{
		PutByte(data_ >> 24);
		PutByte((data_ >> 16) & 0xFF);
		PutByte((data_ >> 8) & 0xFF);
		PutByte(data_ & 0xFF);
	}

	void FlushFullWordBuffer()
	{
		Put4Bytes(wordBuffer);

		wordBuffer    = 0;
		wordBufferPos = 0;
	}

	void FlushPartialWordBuffer()
	{
		wordBuffer <<= (WordBufferSize - wordBufferPos) & 7;

		if (wordBufferPos > 24)
			PutByte(wordBuffer >> 24);
		if (wordBufferPos > 16)
			PutByte((wordBuffer >> 16) & 0xFF);
		if (wordBufferPos > 8)
			PutByte((wordBuffer >> 8) & 0xFF);
		if (wordBufferPos > 0)
			PutByte(wordBuffer & 0xFF);

		wordBuffer = 0;
		wordBufferPos = 0;
	}

	void FillLastByte()
	{
		//os->buf[os->bufPos] <<= (7 - os->bitPos);
		//os->bitPos = 0;
		//os->bufPos += 1;
		uint32 wbp = (WordBufferSize - wordBufferPos) & 7;

		while (wbp > 0)	// can be further optimised
		{
			PutBit(0);
			wbp--;
		}
	}

	void Flush()
	{
		FlushPartialWordBuffer();
	}

	void Reset()
	{
		position = 0;
		wordBufferPos = 0;
		wordBuffer = 0;
	}

private:
	static const uint32 WordBufferSize = 32;

	Buffer&	buffer;
	uint64	position;

	uint32 wordBuffer;
	uint32 wordBufferPos;


	uint32 BitMask(uint32 n_)
	{
		ASSERT(n_ < 32);
		return ((uint32)1 << n_) - 1;
	}
};


#endif // H_BITMEMORY
