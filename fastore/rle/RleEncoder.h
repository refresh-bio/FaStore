/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef RLEENCODER_H
#define RLEENCODER_H

#include "../fastore_bin/Globals.h"
#include "../fastore_bin/BitMemory.h"


/**
 * Binary RLE encoder -- encodes runs of binary values
 *
 */
class BinaryRleEncoder : public ICoder
{
public:
	BinaryRleEncoder(BitMemoryWriter& writer_)
		:	writer(writer_)
		,	currentCount(0)
	{}

	void PutSymbol(bool s_)
	{
		if (s_)
		{
			currentCount++;
			if(currentCount == RleMax - RleOffset)
			{
				Put(currentCount + RleOffset);
				currentCount = 0;
			}
		}
		else
		{
			bool mism = (currentCount > 0) && (currentCount < RleMax - RleOffset);
			if (currentCount > 0)
			{
				Put(currentCount + RleOffset);
				currentCount = 0;
			}
			if (!mism)
				Put(0);
		}
	}

	void Start()
	{
		currentCount = 0;
	}

	void End()
	{
		if (currentCount > 0)
		{
			Put(currentCount + RleOffset);
			currentCount = 0;
		}
	}

private:
	static const uint32 RleMax = 255;
	static const uint32 RleOffset = 2;

	BitMemoryWriter& writer;

	uint32 currentCount;

	void Put(uint32 s_)
	{
		writer.PutByte(s_);
	}
};

class BinaryRleDecoder : public ICoder
{
public:
	BinaryRleDecoder(BitMemoryReader& reader_)
		:	reader(reader_)
		,	currentCount(0)
		,	onlyMatches(false)
	{}

	void Start()
	{
		currentCount = 0;
		onlyMatches = false;
		Fetch();
	}

	void End()
	{}

	bool GetSym()
	{
		if (currentCount > RleOffset)
		{
			currentCount--;
			return true;
		}
		if (currentCount == 0 || (!onlyMatches && currentCount == RleOffset))
		{
			Fetch();
			return false;
		}

		Fetch();
		return GetSym();
	}

private:
	static const uint32 RleOffset = 2;
	static const uint32 RleMax = 255;

	BitMemoryReader& reader;
	uint32 currentCount;
	bool onlyMatches;

	void Fetch()
	{
		if (reader.Position() < reader.Size())
		{
			currentCount = reader.GetByte();
			onlyMatches = (currentCount == RleMax);
		}
	}
};


/**
 * RLE-0 encoder -- encodes 8-, 16- and 32-bit words mixed
 *
 */
class Rle0Encoder : public ICoder
{
public:
	Rle0Encoder(BitMemoryWriter& writer_)
		:	writer(writer_)
		,	prevSymbol(Rle0BSymbol)
	{}

	void PutSymbol(uint32 s_)
	{
		if (s_ == 0)
		{
			if (prevSymbol == Rle0BSymbol)
				prevSymbol = Rle0ASymbol;
			else if (prevSymbol == Rle0ASymbol)
			{
				writer.PutByte(Rle0BSymbol);
				prevSymbol = Rle0BSymbol;
			}
		}
		else
		{
			if (prevSymbol == Rle0ASymbol)
			{
				writer.PutByte(Rle0ASymbol);
				prevSymbol = Rle0BSymbol;
			}

			uint32 ss = s_ + RleOffset;
			if (ss < Max8BitValue)
			{
				writer.PutByte(ss);
			}
			else
			{
				if (ss < (1 << 16) - 1)
				{
					writer.PutByte(Use16Bits);
					writer.Put2Bytes(ss);
				}
				else
				{
					writer.PutByte(Use32Bits);
					writer.Put4Bytes(ss);
				}
			}
		}
	}

	void Start()
	{
		prevSymbol = Rle0BSymbol;
	}

	void End()
	{
		if (prevSymbol == Rle0ASymbol)
			writer.PutByte(Rle0ASymbol);
	}

private:
	static const uint32 Rle0ASymbol = 1;
	static const uint32 Rle0BSymbol = 0;
	static const uint32 RleOffset = 1;

	static const uint32 Max8BitValue = 255 - 2;
	static const uint32 Use16Bits = 0xFE;
	static const uint32 Use32Bits = 0xFF;

	BitMemoryWriter& writer;

	uint32 prevSymbol;
};

class Rle0Decoder : public ICoder
{
public:
	Rle0Decoder(BitMemoryReader& reader_)
		:	reader(reader_)
		,	curSymbol(Rle0BSymbol)
	{}

	void Start()
	{
		curSymbol = Max8BitValue;
	}

	void End()
	{}

	uint32 GetSym()
	{
		if (curSymbol == Rle0BSymbol)
		{
			curSymbol = Rle0ASymbol;
			return 0;
		}
		if (curSymbol == Rle0ASymbol)
		{
			curSymbol = Max8BitValue;
			return 0;
		}

		curSymbol = Fetch();
		if (curSymbol != Rle0ASymbol && curSymbol != Rle0BSymbol)
			return curSymbol - RleOffset;
		else
			return GetSym();
	}

private:
	static const uint32 Rle0ASymbol = 1;
	static const uint32 Rle0BSymbol = 0;
	static const uint32 RleOffset = 1;

	static const uint32 Max8BitValue = 255 - 2;
	static const uint32 Use16Bits = 0xFE;
	static const uint32 Use32Bits = 0xFF;

	BitMemoryReader& reader;
	uint32 curSymbol;

	uint32 Fetch()
	{
		if (reader.Position() < reader.Size())
		{
			uint32 b = reader.GetByte();

			if (b < Use16Bits)
				return b;

			if (b == Use16Bits)
			{
				b = reader.Get2Bytes();
			}
			else
			{
				b = reader.Get4Bytes();
			}

			return b;
		}

		return 0;
	}
};


#endif // RLEENCODER_H
