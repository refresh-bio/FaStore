/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_REBINPARAMS
#define H_REBINPARAMS

#include "../core/Globals.h"

#include <map>

#include "../fastore_bin/Params.h"
#include "../fastore_pack/Params.h"


struct BinBalanceParameters
{
	struct Default : public ReadsClassifierParams::Default, public ReadsContigBuilderParams::Default
	{
		static const uint32 SignatureParity = 2;
		static const uint32 MinBinSizeToExtract = 0;
		static const uint32 MinBinSizeToCategorize = 0;
		static const uint32 MinTreeSize = 4;
		static const bool SelectMaxEdgeRead = true;
		//static const bool FilterReadsWithSignature = false;
	};

	ReadsClassifierParams classifier;
	ReadsContigBuilderParams consensus;

	uint32 signatureParity;
	uint32 minBinSizeToExtract;
	uint32 minBinSizeToCategorize;
	uint32 minTreeSize;
	bool selectMaxEdgeRead;
	//bool filterReadsWithSignature;

	std::vector<bool> validBinSignatures;

	BinBalanceParameters()
		:	signatureParity(Default::SignatureParity)
		,	minBinSizeToExtract(Default::MinBinSizeToExtract)
		,	minBinSizeToCategorize(Default::MinBinSizeToCategorize)
		,	minTreeSize(Default::MinTreeSize)
		,	selectMaxEdgeRead(Default::SelectMaxEdgeRead)
		//,	filterReadsWithSignature(Default::FilterReadsWithSignature)
	{}

	static bool IsSignatureValid(uint32 signatureId_, uint32 signatureParity_)
	{
		// in the context of bins processing:
		// if mod N == 0   then skip		(bins to be merged into - the current iteraton)
		// if mod N/2 == 0 then process		(merged bins from the previous iteration)
		// else skip						(other already processed bins)
		//return ! (  signatureId_ % signatureParity_ == 0 ||
		//		  !(signatureId_ % (signatureParity_ / 2) == 0));
		return  signatureId_ % signatureParity_ != 0 &&
				  signatureId_ % (signatureParity_ / 2) == 0;
	}
};


#endif // H_BINPARAMS
