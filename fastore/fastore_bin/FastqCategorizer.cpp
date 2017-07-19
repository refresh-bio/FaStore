/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"

#include <string.h>
#include <algorithm>

#include "FastqCategorizer.h"


FastqCategorizerBase::FastqCategorizerBase(const MinimizerParameters& params_,
										   const MinimizerFilteringParameters& filter_,
										   const CategorizerParameters& catParams_)
	:	params(params_)
	,	filter(filter_)
	,	catParams(catParams_)
	,	maxLongMinimValue(1 << (2 * params.signatureLen))
	,	nBinValue(maxLongMinimValue)
{
	std::fill(symbolIdxTable.begin(), symbolIdxTable.end(), -1);
	for (uint32 i = 0; i < 5; ++i)
		symbolIdxTable[(int32)params.dnaSymbolOrder[i]] = i;

	InitializeValidBinSignatures(params.signatureLen);
}

void FastqCategorizerBase::InitializeValidBinSignatures(uint32 signatureLen_)
{
	const uint32 totalSignatures = 1 << (signatureLen_ * 2);		
	validBinSignatures.resize(totalSignatures);

	// KMC2 signatures #1
	const uint32 AAA_mask = 0b000000;		// AAA......
	const uint32 AAC_mask = 0b000001;		// AAC......
	const uint32 exc_mask = 0b1111;			// ..AA.....

	// ORCOM signatures
	//const uint32 tri_mask = 0b00111111;
	//const uint32 triplets[] = {0b00000000, 0b00010101, 0b00101010, 0b00111111}; // [A]AAA, [A]CCC, [A]GGG, [A]TTT
	//const uint32 triplets[] = {0x00, 0x15, 0x2A, 0x3F}; // [A]AAA, [A]CCC, [A]GGG, [A]TTT

	const uint32 loMask = (1 << (params.signatureMaskCutoffBits)) - 1;

	for (uint32 i = 0; i < totalSignatures; ++i)
	{
		bool isInvalid = (i & loMask) != 0;

		// KMC2 signatures #2
		uint32 m = i >> (2 * signatureLen_ - 6);
		isInvalid |= (m == AAA_mask) || (m == AAC_mask);

		m = i;
		for (uint32 j = 0; !isInvalid && j < signatureLen_ - 2; ++j)
		{
			// KMC2 signatures
			isInvalid |= ((m & exc_mask) == 0);

			// ORCOM signatures
#if 0
			const uint32 x = m & tri_mask;
			for (uint32 k = 0; k < 4; ++k)
				isInvalid |= (x == triplets[k]);
#endif
			m >>= 2;
		}

		validBinSignatures[i] = !isInvalid;
	}
}


std::pair<uint32, uint16> FastqCategorizerBase::FindMinimizer(const FastqRecord &rec_)
{
	uint32 minimizer = maxLongMinimValue;
	uint16 pos = 0;
	ASSERT(rec_.seqLen >= params.signatureLen - params.skipZoneLen + 1);


	const int32 off0 =  0;
	const int32 off1 = params.skipZoneLen;
	for (int32 i = off0; i < rec_.seqLen - params.signatureLen - off1; ++i)
	{
		// TODO: signature iterator
		uint32 m = ComputeMinimizer(rec_.seq + i, params.signatureLen);

		if (m < minimizer && IsMinimizerValid(m, params.signatureLen))
		{
			minimizer = m;
			pos = i;
		}
	}

	// filter the reads for which we cannot find a proper minimizer
	// and the ones which contain too much N symbols
	if (minimizer >= maxLongMinimValue || std::count(rec_.seq, rec_.seq + rec_.seqLen, 'N') >= rec_.seqLen / 3)
		return std::make_pair(nBinValue, 0);

	return std::make_pair(minimizer, pos);
}


std::map<uint32, uint16> FastqCategorizerBase::FindMinimizers(const FastqRecord &rec_, uint32 startOff_, uint32 endCutoff_)
{
	ASSERT(rec_.seqLen >= params.signatureLen - params.skipZoneLen + 1);

	const int32 off0 = startOff_;
	const int32 off1 = params.skipZoneLen + endCutoff_;

	// find all signatures
	//
	std::map<uint32, uint16> signatures;
	for (int32 i = off0; i < rec_.seqLen - params.signatureLen - off1; ++i)
	{
		// TODO: signature iterator
		uint32 m = ComputeMinimizer(rec_.seq + i, params.signatureLen);

		if (m < maxLongMinimValue
				&& IsMinimizerValid(m, params.signatureLen)
				&& (!filter.filterLowQualitySignatures
					|| (rec_.qua != NULL && IsMinimizerQualityValid(rec_.qua + i, params.signatureLen, filter.lowQualityThreshold))))
		{
			if (signatures.count(m) == 0)
				signatures[m] = i + (params.signatureLen - params.signatureLen);
		}
	}

	return signatures;
}


uint32 FastqCategorizerBase::ComputeMinimizer(const char* dna_, uint32 mLen_)
{
	uint32 r = 0;

	for (uint32 i = 0; i < mLen_; ++i)
	{
		if (dna_[i] == 'N')
			return nBinValue;

		ASSERT(dna_[i] >= 'A' && dna_[i] <= 'T');
		r = (r << 2) + symbolIdxTable[(uint32)dna_[i]];
	}

	return r;
}

bool FastqCategorizerBase::IsMinimizerValid(uint32 minim_, uint32 /*mLen_*/)
{
	ASSERT(minim_ != maxLongMinimValue);
	return validBinSignatures[minim_];
}

bool FastqCategorizerBase::IsMinimizerQualityValid(const char *qua_, uint32 mLen_, uint32 minQ_)
{
	bool valid = true;
	for (uint32 i = 0; i < mLen_; ++i)
		valid = valid && ((uint32)qua_[i] >= minQ_);

	return valid;
}

void FastqCategorizerSE::Categorize(std::vector<FastqRecord>& records_,
									std::map<uint32, FastqRecordsPtrBin>& bins_)
{
	ASSERT(!records_.empty());

	// crear bins
	//
	bins_.clear();

	// process records
	//
	DistributeToBins(records_, bins_);

	// sort the bins
	//
	/*
	TFastqComparator<const FastqRecordSE*> comparator;
	for (FastqPtrBinMapSE::iterator i = bins_.begin(); i != bins_.end(); ++i)
	{
		std::vector<FastqRecordSE*>& db = i->second.records;
		std::sort(db.begin(), db.end(), comparator);
	}
	*/
}


// TODO: split into rev and non-rev
//
void FastqCategorizerSE::DistributeToBins(std::vector<FastqRecord>& records_,
										  std::map<uint32, FastqRecordsPtrBin>& bins_)
{
	FastqRecordBuffer rcRec;

	for (FastqRecord& rec : records_)
	{
		rec.SetReadReverse(false);

		ASSERT(rec.seqLen > 0);

		rec.ComputeRC(rcRec);

		// find and select minimizers
		//
		auto minimizerFwd = FindMinimizer(rec);
		auto minimizerRev = FindMinimizer(rcRec);
		decltype(minimizerFwd) minimizer;
		bool reverse = false;

		if (minimizerFwd.first <= minimizerRev.first)
		{
			minimizer = minimizerFwd;
		}
		else
		{
			minimizer = minimizerRev;
			reverse = true;
		}

		// store record to bin
		//
		FastqRecordsPtrBin* rb = NULL;
		if (minimizer.first != nBinValue)								// !TODO --- find here minimizer pos
		{
			rb = &bins_[minimizer.first];
			if (reverse)
			{
				rec.SetReadReverse(true);
                                rec.CopyFrom(rcRec);
			}
			rec.minimPos = minimizer.second;
		}
		else
		{
			rb = &bins_[nBinValue];
			rec.minimPos = 0;
			rec.SetReadReverse(false);
		}

		rb->records.push_back(&rec);
		rb->stats.Update(rec);

		ASSERT(rb->stats.minSeqLen > 0);
		ASSERT(rb->stats.maxSeqLen);
	}
}


void FastqCategorizerPE::DistributeToBins(std::vector<FastqRecord> &records_,
										  std::map<uint32, FastqRecordsPtrBin> &bins_)
{
	FastqRecordBuffer recRev;

	for (FastqRecord& rec : records_)
	{
		rec.Reset();

		ASSERT(rec.seqLen > 0);
		ASSERT(rec.auxLen > 0);

		recRev.seqLen = rec.seqLen;
		recRev.auxLen = rec.auxLen;

		const FastqRecord rec_2 = rec.GetPair();
		rec.ComputeRC(recRev);

		const FastqRecord recRev_2 = recRev.GetPair();

		// find and select minimizers
		//
		bool isRev = false;
		bool isFwdMinim = true;

		auto minFwd_1 = FindMinimizer(rec);
		auto minRev_1 = FindMinimizer(recRev);			// old _2
		decltype(minFwd_1) minimizer;

#if 1
		auto minFwd_2 = FindMinimizer(rec_2);
		auto minRev_2 = FindMinimizer(recRev_2);

		bool isFwdMinim_1 = minFwd_1.first < minFwd_2.first;
		auto fwdMinim = isFwdMinim_1 ? minFwd_1 : minFwd_2;

		bool isRevMinim_1 = minRev_1.first < minRev_2.first;
		auto revMinim = isRevMinim_1 ? minRev_1 : minRev_2;

		if (fwdMinim.first < revMinim.first)
		{
			minimizer = fwdMinim;
			isFwdMinim = isFwdMinim_1;
		}
		else
		{
			minimizer = revMinim;
			isRev = true;
			isFwdMinim = isRevMinim_1;
		}
#else
		if (minFwd_1.first < minRev_1.first)
		{
			minimizer = minFwd_1;
		}
		else
		{
			minimizer = minRev_1;
			isRev = true;
		}
#endif

		// store record to bin
		//
		FastqRecordsPtrBin* rb = NULL;
		if (minimizer.first != nBinValue)
		{
			rb = &bins_[minimizer.first];
			if (isRev)
			{
				rec.CopyFrom(recRev);
				rec.SetReadReverse(true);
			}

			// we swap here the reads R1 <-> R2
			//
			// TODO: we can optimize here to avoid extra copy in case of rev-compl and swapping
			if (!isFwdMinim)
				rec.SwapReads();

			rec.minimPos = minimizer.second;

#if DEV_DEBUG_MODE
		char minString[64];
		params.GenerateMinimizer(minimizer.first, minString);

		const char* minSeq = std::search(rec.seq,
										 rec.seq + rec.seqLen,
										 minString,
										 minString + params.signatureLen);

		ASSERT(minSeq != rec.seq + rec.seqLen);
		ASSERT(minSeq == rec.seq + rec.minimPos);
#endif

		}
		else
		{
			rb = &bins_[nBinValue];
		}

		rb->records.push_back(&rec);
		rb->stats.maxSeqLen = MAX(rb->stats.maxSeqLen, rec.seqLen);
		rb->stats.minSeqLen = MIN(rb->stats.minSeqLen, rec.seqLen);
		rb->stats.maxAuxLen = MAX(rb->stats.maxAuxLen, rec.auxLen);
		rb->stats.minAuxLen = MIN(rb->stats.minAuxLen, rec.auxLen);
	}
}
