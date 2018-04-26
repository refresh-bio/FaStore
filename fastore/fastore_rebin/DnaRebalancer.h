/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DNAREBALANCER
#define H_DNAREBALANCER

#include "../core/Globals.h"

#include <vector>
#include <tuple>

#include "Params.h"

#include "../core/Node.h"
#include "../core/FastqCategorizer.h"
#include "../core/ReadsClassifier.h"


/**
 * Re-bins the binned reads
 *
 */
class DnaRebalancer : public FastqCategorizerBase
{
public:
	DnaRebalancer(const MinimizerParameters& params_,
				  const BinBalanceParameters& binParams_,
				  bool pairedEnd_ = false);

	void Rebalance(RebinContext& rebinCtx_,
				   std::map<uint32, MatchNodesPtrBin>& bins_,
				   uint32 signature_);

protected:
	const BinBalanceParameters binParams;
	const bool pairedEnd;

	ReadsClassifierSE readsClassifier;

	void StoreSingleReadNode(MatchNode* node_,
							 std::map<uint32, MatchNodesPtrBin>& bins_,
							 uint32 signature_);

	void StoreTree(MatchNode* root_,
				   RebinContext& ctx_,
				   std::map<uint32, MatchNodesPtrBin>& bins_,
				   uint32 signature_);

	std::pair<const MatchNode*, const MatchNode*> FindEdgeNodes(const MatchNode* root_);


	using FastqCategorizerBase::FindMinimizer;
	using FastqCategorizerBase::FindMinimizers;

	std::pair<uint32, uint16> FindMinimizerHR(const FastqRecord &rec_,
											  uint32 curSignature_,
											  uint32 curDivisor_);

	std::tuple<uint32, uint16, bool> FindNewMinimizer(const FastqRecord &recFwd_,
													  const FastqRecord &recRev_,
													  uint32 curSignature_);

	void UpdateExactMatches(MatchNode* node_, bool directionChange_);
	void UpdateTreeReads(MatchNode* node_, bool directionChange_);
};


#endif // H_DNAREBALANCER
