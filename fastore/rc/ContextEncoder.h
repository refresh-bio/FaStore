/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_CONTEXTENCODER
#define H_CONTEXTENCODER

#include "../fastore_bin/Globals.h"

#include "RangeCoder.h"
#include "SymbolCoderRC.h"


template <uint32 _TSymbolCount, uint32 _TSymbolOrder, uint32 _TTotalOrder = _TSymbolOrder>
class TDynamicContextCoderBase
{
public:
	static const uint32 SymbolCount = _TSymbolCount;
	static const uint32 SymbolOrder = _TSymbolOrder;
	static const uint32 TotalOrder = _TTotalOrder;

	static const uint32 AlphabetSize = SymbolCount;

	TDynamicContextCoderBase()
		:	models(NULL)
		,	hash(0)
		,	symBuffer(0)
	{
		models = new Coder[ModelCount];
	}

	~TDynamicContextCoderBase()
	{
		delete models;
	}

	void Clear()
	{
		hash = 0;
		symBuffer = 0;
		CoderStatType* cd = (CoderStatType*)models;
		std::fill(cd, cd + ModelCount * SymbolCount, 1);
	}

protected:
	typedef uint64 THash;

	static const uint32 ModelCount = TPow<AlphabetSize, TotalOrder>::Value;

	static const uint64 HiPow = TPow<AlphabetSize, _TSymbolOrder-1>::Value;

	typedef TSymbolCoderRC<AlphabetSize> Coder;
	typedef typename Coder::StatType CoderStatType;

	Coder* models;
	THash hash;
	THash symBuffer;

	void UpdateHash(uint32 sym_)
	{
		if (ModelCount == 1)
			return;

		hash = (hash - HiPow * (hash / HiPow)) * AlphabetSize + sym_;
		ASSERT(hash < ModelCount);
	}

	THash GetHash() const
	{
		if (ModelCount == 1)
			return 1;

		return hash;
	}
};


template <uint32 _TAlphabetSize, uint32 _TSymbolOrder, uint32 _TTotalOrder = _TSymbolOrder>
class TStaticContextCoderBase
{
public:
	static const uint32 AlphabetSize = _TAlphabetSize;
	static const uint32 AlphabetBits = TLog2<AlphabetSize>::Value;
	static const uint32 SymbolOrder = _TSymbolOrder;
	static const uint32 TotalOrder = _TTotalOrder;

	TStaticContextCoderBase()
		:	hash(0)
	{}

	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_)
	{
		models[GetHash()].EncodeSymbol(rc_, sym_);

		UpdateHash(sym_);
	}

	uint32 DecodeSymbol(RangeDecoder& rc_)
	{
		uint32 sym = models[GetHash()].DecodeSymbol(rc_);

		UpdateHash(sym);

		return sym;
	}

	void Clear()
	{
		// clear hash
		hash = 0;

		// clear stats -- fill context data with initial '1' value
		CoderStatType* cd = (CoderStatType*)models;
		std::fill(cd, cd + ModelCount * AlphabetSize, 1);
	}

protected:
	typedef uint64 HashType;
	typedef TSymbolCoderRC<AlphabetSize> Coder;
	typedef typename Coder::StatType CoderStatType;

	static const HashType HashMask = (1ULL << (TotalOrder * AlphabetBits)) - 1;
	static const HashType SymbolHashMask = (1ULL << (SymbolOrder * AlphabetBits)) - 1;

	static const uint32 ModelCount = 1 << (TLog2<AlphabetSize>::Value * TotalOrder);

	Coder models[ModelCount];
	HashType hash;

	HashType GetHash()
	{
		return hash & SymbolHashMask;
	}

	void UpdateHash(uint32 sym_)
	{
		hash <<= AlphabetBits;
		hash |= sym_;
		//hash &= HashMask;
	}
};

template <uint32 _TSymbolCount, uint32 _TOrder>
class TSimpleContextCoder : public TStaticContextCoderBase<_TSymbolCount, _TOrder>
{
public:
	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_)
	{
		Super::models[Super::GetHash()].EncodeSymbol(rc_, sym_);

		Super::UpdateHash(sym_);
	}

	uint32 DecodeSymbol(RangeDecoder& rc_)
	{
		uint32 sym = Super::models[Super::GetHash()].DecodeSymbol(rc_);

		Super::UpdateHash(sym);
		return sym;
	}

private:
	typedef TStaticContextCoderBase<_TSymbolCount, _TOrder> Super;
};


template <uint32 _TSymbolCount, uint32 _TOrder>
class TAdvancedContextCoder : public TStaticContextCoderBase<_TSymbolCount, _TOrder, _TOrder + 1>
{
public:
	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_, uint32 ctx0_ = 0)
	{
		ASSERT(ctx0_ < Super::AlphabetSize);
		ASSERT(sym_ < Super::AlphabetSize);

		uint32 h = (Super::GetHash() << Super::AlphabetBits) | ctx0_;
		//uint64 h = Super::GetHash() + ctx0_ * HiCtxPow;
		ASSERT(h < Super::ModelCount);

		Super::models[h].EncodeSymbol(rc_, sym_);

		Super::UpdateHash(sym_);
	}

	uint32 DecodeSymbol(RangeDecoder& rc_, uint32 ctx0_ = 0)
	{
		ASSERT(ctx0_ < Super::AlphabetSize);

		uint32 h = (Super::GetHash()  << Super::AlphabetBits) | ctx0_;
		//uint64 h = Super::GetHash() + ctx0_ * HiCtxPow;
		ASSERT(h < Super::ModelCount);

		uint32 sym = Super::models[h].DecodeSymbol(rc_);

		Super::UpdateHash(sym);
		return sym;
	}

private:
	typedef TStaticContextCoderBase<_TSymbolCount, _TOrder, _TOrder + 1> Super;
};


template <class _TRangeCoder, class _TBitMemory>
struct TCoderBase : public ICoder
{
	_TRangeCoder rc;

	TCoderBase(_TBitMemory& mem_)
		:	rc(mem_)
	{}

	void Start()
	{
		rc.Start();
	}

	void End()
	{
		rc.End();
	}
};


template <class _TContextEncoder>
struct TEncoder : public TCoderBase<RangeEncoder, BitMemoryWriter>
{
	typedef TCoderBase<RangeEncoder, BitMemoryWriter> Super;

	_TContextEncoder coder;

	TEncoder(BitMemoryWriter& mem_)
		: Super(mem_)
	{}

	void Start()
	{
		Super::Start();
		coder.Clear();
	}
};


template <class _TContextEncoder>
struct TDecoder : public TCoderBase<RangeDecoder, BitMemoryReader>
{
	typedef TCoderBase<RangeDecoder, BitMemoryReader> Super;

	_TContextEncoder coder;

	TDecoder(BitMemoryReader& mem_)
		: Super(mem_)
	{}

	void Start()
	{
		Super::Start();
		coder.Clear();
	}
};


#endif // H_CONTEXTENCODER
