/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"
#include "FastqStream.h"
#include "Utils.h"


uint64 IFastqStreamReaderBase::GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_)
{
	SkipToEol(data_, pos_, size_);
	++pos_;

	// find beginning of the next record
	//
	while (data_[pos_] != '@')
	{
		SkipToEol(data_, pos_, size_);
		++pos_;
	}
	uint64 pos0 = pos_;

	SkipToEol(data_, pos_, size_);
	++pos_;

	if (data_[pos_] == '@')				// previous one was a quality field
		return pos_;

	SkipToEol(data_, pos_, size_);
	++pos_;

	ASSERT(data_[pos_] == '+');			// pos0 was the start of tag
	return pos0;
}



bool IFastqStreamReaderSE::ReadNextChunk(IFastqChunkCollection& chunk_)
{
	ASSERT(chunk_.chunks.size() == 1);

	if (Eof())
	{
		chunk_.chunks[0]->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk
	//
	uchar* data = chunk_.chunks[0]->data.Pointer();
	const uint64 cbufSize = chunk_.chunks[0]->data.Size();
	chunk_.chunks[0]->size = 0;

	int64 toRead = cbufSize - readBufferSize;
	if (readBufferSize > 0)
	{
		std::copy(readBuffer.Pointer(), readBuffer.Pointer() + readBufferSize, data);
		chunk_.chunks[0]->size = readBufferSize;
		readBufferSize = 0;
	}

	// read the next chunk
	//
	int64 r = Read(data + chunk_.chunks[0]->size, toRead);
	if (r > 0)
	{
		if (r == toRead)				// somewhere before end
		{
			uint64 chunkEnd = cbufSize - maxReadBufferSize;

			chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);

			chunk_.chunks[0]->size = chunkEnd - 1;
			if (usesCrlf)
				chunk_.chunks[0]->size -= 1;

			std::copy(data + chunkEnd, data + cbufSize, readBuffer.Pointer());
			readBufferSize = cbufSize - chunkEnd;
		}
		else							// at the end of file
		{
			chunk_.chunks[0]->size += r - 1;		// skip the last EOF symbol
			if (usesCrlf)
				chunk_.chunks[0]->size -= 1;

			eof = true;
		}
	}
	else
	{
		eof = true;
	}

	return true;
}


bool IFastqStreamReaderPE::ReadNextChunk(IFastqChunkCollection& chunk_)
{
	ASSERT(chunk_.chunks.size() >= 2);

	if (Eof())
	{
		chunk_.chunks[0]->size = 0;
		chunk_.chunks[1]->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk : _1
	//
	uchar* data_1 = chunk_.chunks[0]->data.Pointer();
	const uint64 cbufsz_1 = chunk_.chunks[0]->data.Size();
	chunk_.chunks[0]->size = 0;

	int64 toRead_1 = cbufsz_1 - readBufferSize;
	if (readBufferSize > 0)
	{
		std::copy(readBuffer.Pointer(), readBuffer.Pointer() + readBufferSize, data_1);
		chunk_.chunks[0]->size = readBufferSize;
		readBufferSize = 0;
	}

	// flush the data from previous incomplete chunk : _2
	//
	uchar* data_2 = chunk_.chunks[1]->data.Pointer();
	const uint64 cbufsz_2 = chunk_.chunks[1]->data.Size();
	chunk_.chunks[1]->size = 0;

	int64 toRead_2 = cbufsz_2 - pairBufferSize;
	if (pairBufferSize > 0)
	{
		std::copy(pairBuffer.Pointer(), pairBuffer.Pointer() + pairBufferSize, data_2);
		chunk_.chunks[1]->size = pairBufferSize;
		pairBufferSize = 0;
	}

	// read the next chunks : _1 and _2
	//
	int64 r_1 = Read_1(data_1 + chunk_.chunks[0]->size, toRead_1);
	int64 r_2 = Read_2(data_2 + chunk_.chunks[1]->size, toRead_2);

	// synchronise the chunks
	//
	if (r_1 == toRead_1 && r_2 == toRead_2)
	{
		ASSERT(cbufsz_1 >= maxReadBufferSize);
		ASSERT(cbufsz_2 >= maxReadBufferSize);

		uint64 chunkEnd_1 = cbufsz_1 - maxReadBufferSize;
		uint64 chunkEnd_2 = cbufsz_2 - maxReadBufferSize;

		chunkEnd_1 = GetNextRecordPos(data_1, chunkEnd_1, cbufsz_1);
		chunkEnd_2 = GetNextRecordPos(data_2, chunkEnd_2, cbufsz_2);

		// get the read IDs
		//
		uint64 rid_1 = ParseNextReadId(data_1 + chunkEnd_1, maxReadBufferSize);
		uint64 rid_2 = ParseNextReadId(data_2 + chunkEnd_2, maxReadBufferSize);

		// those ones should be usually the same, synchronised
		//

		if (rid_1 < rid_2)		// read more records _1
		{
			const uint64 nl = (rid_2 - rid_1) * 4;
			for (uint64 i = 0; i < nl; ++i)
			{
				SkipToEol(data_1, chunkEnd_1, cbufsz_1);
				chunkEnd_1++;
				rid_1++;
			}
		}
		if (rid_1 > rid_2)		// read more records _2
		{
			const uint64 nl = (rid_1 - rid_2) * 4;
			for (uint64 i = 0; i < nl; ++i)
			{
				SkipToEol(data_2, chunkEnd_2, cbufsz_2);
				chunkEnd_2++;
				rid_2++;
			}
		}

		ASSERT(rid_1 == rid_2);

		chunk_.chunks[0]->size = chunkEnd_1 - 1;
		if (usesCrlf)
			chunk_.chunks[0]->size -= 1;

		std::copy(data_1 + chunkEnd_1, data_1 + cbufsz_1, readBuffer.Pointer());
		readBufferSize = cbufsz_1 - chunkEnd_1;

		chunk_.chunks[1]->size = chunkEnd_2 - 1;
		if (usesCrlf)
			chunk_.chunks[1]->size -= 1;

		std::copy(data_2 + chunkEnd_2, data_2 + cbufsz_2, pairBuffer.Pointer());
		pairBufferSize = cbufsz_2 - chunkEnd_2;
	}
	else
	{
		ASSERT((r_1 > 0 && r_2 > 0) || (r_1 == 0 && r_2 == 0));

		if (r_1 > 0)
		{
			chunk_.chunks[0]->size += r_1 - 1;		// skip the last EOF symbol
			if (usesCrlf)
				chunk_.chunks[0]->size -= 1;
		}
		eof = true;

		if (r_2 > 0)
		{
			chunk_.chunks[1]->size += r_2 - 1;		// skip the last EOF symbol
			if (usesCrlf)
				chunk_.chunks[1]->size -= 1;
		}
		eof_2 = true;
	}

	return true;
}


uint64 IFastqStreamReaderPE::ParseNextReadId(const uchar *data_, uint64 maxLen_)
{
	const char *sep = " ._,=:/-#"; //9
	const std::vector<uchar> separators(sep, sep + 9 + 1);

	const uchar* tag = NULL;
	uint64 len = 0;

	for (const uchar* p = data_; p < data_ + maxLen_; ++p)
	{
		if (!std::count(separators.begin(), separators.end(), *p) && (*p != '\n'))
			continue;

		if (tag == NULL)
		{
			tag = ++p;
		}
		else
		{
			len = (uint64)(p - tag);
			break;
		}
	}

	return to_num(tag, len);
}
