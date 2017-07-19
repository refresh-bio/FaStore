/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_SYMBOLCODERRC
#define H_SYMBOLCODERRC

#include "../fastore_bin/Globals.h"

#include "RangeCoder.h"


template <uint32 _TMaxSymbolCount>
class TSymbolCoderRC
{
public:
	typedef uint16 StatType;
	static const uint32 MaxSymbolCount = _TMaxSymbolCount > 0 ? _TMaxSymbolCount : 1;

	TSymbolCoderRC()
	{
		Clear();
	}

	void EncodeSymbol(RangeEncoder& rc_, uint32 sym_)
	{
		ASSERT(sym_ < MaxSymbolCount);

		uint32 acc = Accumulate();
		uint32 loEnd = 0;

		for (uint32 i = 0; i < sym_; ++i)
			loEnd += stats[i];

		rc_.EncodeFrequency(stats[sym_], loEnd, acc);

		stats[sym_] += StepSize;
	}

	uint32 DecodeSymbol(RangeDecoder& rc_)
	{
		uint32 acc = Accumulate();
		uint32 cul = rc_.GetCumulativeFreq(acc);

		uint32 idx, hiEnd;
		for (idx = 0, hiEnd = 0; (hiEnd += stats[idx]) <= cul; ++idx)
			;
		hiEnd -= stats[idx];

		rc_.UpdateFrequency(stats[idx], hiEnd, acc);
		stats[idx] += StepSize;
		return idx;
	}

	void Clear()
	{
		std::fill(stats, stats + MaxSymbolCount, 1);
	}

private:
	static const StatType StepSize = 8;
	static const uint32 MaxAccumulatedValue = (1<<16) - MaxSymbolCount*StepSize;

	void Rescale()
	{
		for (uint32 i = 0; i < MaxSymbolCount; ++i)
			stats[i] -= stats[i] >> 1;		// no '>>=' to avoid reducing stats to 0
	}

	uint32 Accumulate()
	{
		uint32 acc = 0;
		for (uint32 i = 0; i < MaxSymbolCount; ++i)
			acc += stats[i];

		if (acc >= MaxAccumulatedValue)
		{
			acc = 0;
			Rescale();
			for (uint32 i = 0; i < MaxSymbolCount; ++i)
				acc += stats[i];
		}

		return acc;
	}

	StatType stats[MaxSymbolCount];			// can be assumed to be uint16
};


#endif // H_SYMBOLCODERRC
