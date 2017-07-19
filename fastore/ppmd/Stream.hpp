/****************************************************************************
 *  This file is part of PPMd library                                       *
 *  Contents: Memory Byte Stream wrapper for I/O operations                 *
 *  Comments: exchanges the raw I/O routines used in original PPMd project  *
 *  Author: Lucas Roguski, 2014                                             *
 ****************************************************************************/

#ifndef _BYTESTREAM_H_
#define _BYTESTREAM_H_

#include <stdint.h>

#define USE_64BIT_SIZE 1

class ByteStream
{
public:
	ByteStream(unsigned char* memory_, uint64_t size_)
		:	memory(memory_)
		,	size(size_)
		,	position(0)
	{}

	int Get()
	{
		if (position >= size)
			return -1;
		return memory[position++];
	}

	int Put(int c_)
	{
		if (position >= size)
			return -1;
		memory[position++] = (unsigned char)c_;
		return c_;
	}

	unsigned char* Memory()
	{
		return memory;
	}

	uint64_t Size() const
	{
		return size;
	}

	uint64_t Position() const
	{
		return position;
	}

private:
	unsigned char* 	memory;
	uint64_t		size;
	uint64_t		position;
};

#endif // _BYTESTREAM_H_
