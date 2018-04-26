/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef CONSENSUSBUILDER_H
#define CONSENSUSBUILDER_H

#include "../core/Globals.h"

#include <vector>
#include <map>
#include <algorithm>

#include "../core/FastqRecord.h"
#include "../core/ReadsClassifier.h"


/**
 * Auxiliary contig information used when assembling the reads
 *
 */
struct ContigBuildInfo
{
	struct WorkNode
	{
		MatchNode* match;
		std::vector<uint16> newVariantPositions;

		WorkNode()
			:	match(NULL)
		{}

		bool operator<(const WorkNode& r2_) const
		{
			return match->record->minimPos < r2_.match->record->minimPos;
		}
	};

	static const char EmptyChar = '.';
	static const uint32 DefaultSequenceLen = 100;

	ConsensusDefinition consensus;

	// temporary build info
	MatchNode* mainNode;
	std::vector<WorkNode> nodes;

	std::vector<uint16> variantFreqPerPos;
	std::vector<uint16> recordsPerPos;

	std::vector<MatchNode*> removedNodes;
	std::vector<MatchNode*> unlinkedParents;		// DEBUG ONLY


	ContigBuildInfo(uint32 seqLen_ = DefaultSequenceLen)
	{
		ASSERT(seqLen_ != 0);

		Reset(seqLen_);
	}

	void Reset(uint32 seqLen_ = DefaultSequenceLen)
	{
		if (consensus.readLen != seqLen_ || consensus.sequence.size() != seqLen_ * 2)
		{
			consensus.sequence.resize(seqLen_ * 2);
			consensus.variantPositions.resize(seqLen_ * 2);
			consensus.readLen = seqLen_;

			variantFreqPerPos.resize(seqLen_ * 2);
			recordsPerPos.resize(seqLen_ * 2);
		}

		std::fill(consensus.sequence.begin(), consensus.sequence.end(), +EmptyChar);
		std::fill(consensus.variantPositions.begin(), consensus.variantPositions.end(), 0);

		std::fill(variantFreqPerPos.begin(), variantFreqPerPos.end(), 0);
		std::fill(recordsPerPos.begin(), recordsPerPos.end(), 0);

		mainNode = NULL;
		nodes.clear();

		removedNodes.clear();
		unlinkedParents.clear();	

		consensus.variantsCount = 0;
		consensus.range.first = seqLen_;
		consensus.range.second = seqLen_;
	}
};


struct ReadsContigBuilderParams
{
	struct Default
	{
		static const uint32 BeginCut = 2;
		static const uint32 EndCut = 2;
		static const uint32 MaxNewVariantsPerRead = 1;
		static const uint32 MaxRecordShiftDiff = 0;			// INFO: 0 - auto, half of the length (tests: 42)
		static const uint32 MaxHammingDistance = 8;			// TODO: a better metric?
		static const uint32 MinConsensusSize = 10;
	};

	uint32 beginCut;
	uint32 endCut;
	uint32 maxNewVariantsPerRead;
	uint32 maxRecordShiftDifference;
	uint32 maxHammingDistance;
	uint32 minConsensusSize;

	ReadsContigBuilderParams()
		:	beginCut(Default::BeginCut)
		,	endCut(Default::EndCut)
		,	maxNewVariantsPerRead(Default::MaxNewVariantsPerRead)
		,	maxRecordShiftDifference(Default::MaxRecordShiftDiff)
		,	maxHammingDistance(Default::MaxHammingDistance)
		,	minConsensusSize(Default::MinConsensusSize)
	{}
};


class ContigBuilder
{
public:
	ContigBuilder(const ReadsContigBuilderParams& params_,
					 const MinimizerParameters& minParams_);
	~ContigBuilder();

	// returns the (possibly modified) tree root
	bool Build(MatchNode* root_, PackContext& packCtx_);

private:
	static const bool AvoidTreesInConsensus = true;

	const ReadsContigBuilderParams params;
	const MinimizerParameters minParams;


	bool AddRecord(ContigBuildInfo& cons_, MatchNode* node_, bool tryFullMatchOnly_ = false);

	// returns the number of removed records
	uint32 OptimizeContig(ContigBuildInfo& cons_);

	void UpdateContigLinkage(ContigBuildInfo& contig_);

	void PostProcessContig(ContigBuildInfo& cons_);

	void StoreContig(ContigBuildInfo& contig_, PackContext& packCtx_);



	// TODO: explore and tweak better those methods
	//
	float GetNormalEncodeCost(const MatchNode& node_)
	{
		// (mismatch + shift encode cost) + rle cost + lz_id cost
		float rleCost = 0.0f;
		if (ABS(node_.shiftValue) != node_.encodeCost)
			rleCost = 1.0f + node_.encodeCost / 1.5f;

		return (1.0f + node_.encodeCost) + rleCost + 2.0f;
	}

	float GetConsEncodeCost(const ContigBuildInfo& cons_, const ContigBuildInfo::WorkNode& node_, uint32 hamDistance_)
	{
		uint32 newVarCost = 0;
		for (std::vector<uint16>::const_iterator i = node_.newVariantPositions.begin();
			 i != node_.newVariantPositions.end(); ++i)
		{
			newVarCost += cons_.recordsPerPos[*i];
		}

		if (newVarCost > 0)
			hamDistance_ -= 1;

		return (float)(1 + hamDistance_ + params.beginCut + params.endCut) + float(newVarCost) * 0.9;
	}
};


#endif // CONSENSUSBUILDER_H
