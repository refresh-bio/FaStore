/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef READSCLASSIFIER_H
#define READSCLASSIFIER_H

#include "../core/FastqRecord.h"

#include <vector>
#include <deque>
#include <stack>
#include <list>
#include <type_traits>

#include "Node.h"
#include "FastqCategorizer.h"


template <typename _TNode, bool _TConst = false>
struct TNodeIterator
{
	typedef typename std::conditional<_TConst, const _TNode, _TNode>::type Node;
	std::deque<Node> nodes;

	TNodeIterator(Node main_ = NULL)
	{
		nodes.push_back(main_);
	}

	bool Empty() const
	{
		return nodes.empty();
	}

	void Skip()
	{
		ASSERT(!nodes.empty());

		Node n = nodes.front();
		nodes.pop_front();

		if (n->HasChildren())
		{
			EnqueueChildren(n);
		}
	}

	Node Next()
	{
		ASSERT(!nodes.empty());

		Node n = nodes.front();
		nodes.pop_front();

		if (n->HasChildren())
		{
			EnqueueChildren(n);
		}

		return n;
	}

	void EnqueueChildren(Node node_)
	{
		ASSERT(node_->children != NULL);

#if 0
		if (node_->children->size() > 1)
		{
			std::vector<Node> cnodes;
			for (auto child : *node_->children)
			{
				if (child->HasChildren())
					cnodes.push_back(child);
				else
					nodes.push_back(child);
			}
			nodes.insert(nodes.end(), cnodes.begin(), cnodes.end());
		}
		else
		{
			nodes.push_back(*node_->children->begin());
		}
#else
		nodes.insert(nodes.end(), node_->children->begin(), node_->children->end());
#endif
	}
};

typedef TNodeIterator<MatchNode*, false> MatchNodeIterator;
typedef TNodeIterator<MatchNode*, true> MatchCNodeIterator;



struct ReadsClassifierParams
{
	struct Default
	{
		static const int32 MaxCostValue = (uint16)-1;
		static const int32 AutoEncodeThresholdValue = 0;
		static const int32 ShiftCost = 1;
		static const int32 MismatchCost = 2;
		static const uint32 MaxLzWindowSize = MAX_LZ_SE;
		static const uint32 MaxPairLzWindowSize = MAX_LZ_PE;
		static const bool ExtraReduceHardReads = false;			// temporary, for backwards-compatibility
		static const bool ExtraReduceExpensiveLzMatches = false;
	};

	int32 maxCostValue;
	int32 encodeThresholdValue;
	int32 pairEncodeThresholdValue;
	int32 shiftCost;
	int32 mismatchCost;
	uint32 maxLzWindowSize;
	uint32 maxPairLzWindowSize;
	bool extraReduceHardReads;
	bool extraReduceExpensiveLzMatches;

	ReadsClassifierParams()
		:	maxCostValue(Default::MaxCostValue)
		,	encodeThresholdValue(Default::AutoEncodeThresholdValue)
		,	pairEncodeThresholdValue(Default::AutoEncodeThresholdValue)
		,	shiftCost(Default::ShiftCost)
		,	mismatchCost(Default::MismatchCost)
		,	maxLzWindowSize(Default::MaxLzWindowSize)
		,	maxPairLzWindowSize(Default::MaxPairLzWindowSize)
		,	extraReduceHardReads(Default::ExtraReduceHardReads)
		,	extraReduceExpensiveLzMatches(Default::ExtraReduceExpensiveLzMatches)
	{}
};


/**
 * Reads classifers used when performing matching
 *
 */
class ReadsClassifierSE
{
public:
	struct MatchCost
	{
		int32 cost;
		bool noMismatches;

		MatchCost()
			:	cost(255)
			,	noMismatches(false)
		{}
	};

	struct MatchResult
	{
		static const uint32 MaxCost = 255;
		static const uint32 MaxInsert = 128 - 1;

		MatchCost cost;
		int32 prevId;
		int32 shift;

		MatchResult()
			:	prevId(0)
			,	shift(0)
		{}
	};

	// TODO: move out into external LZ matcher
	struct LzMatch
	{
                char* seq;
                MatchNode* node;
                uint16 seqLen;
		uint16 minPos;

		LzMatch()
                        :	seq(NULL)
                        ,	node(NULL)
                        ,	seqLen(0)
			,	minPos(0)
		{}
	};

	ReadsClassifierSE(const MinimizerParameters& minParams_,
					  const ReadsClassifierParams& classifierParams_ = ReadsClassifierParams());
	~ReadsClassifierSE();

	void ConstructMatchTree(GraphEncodingContext& graph_,
							std::vector<MatchNode*>& outRootNodes_,
							MatchNode* auxRootNode_ = NULL);


	// trying to inline this guy for optimization....
	bool UpdateLzMatchResult(MatchResult& result_,
                                 const char* seq_, uint32 seqLen_, int32 minPos_,
								 const char* lzSeq_, uint32 lzSeqLen_, int32 lzMinPos_,
								 uint32 minCostThreshold_ = 0) const
	{
		(void)minCostThreshold_;

		const int32 shift = lzMinPos_ - minPos_;
		const int32 insertCost = ABS(shift) * classifierParams.shiftCost;

		if (insertCost > result_.cost.cost || (uint32)ABS(shift) > MatchResult::MaxInsert)
			return false;

		const int32 recOff = (shift < 0) ? (-shift) : 0;
		const int32 lzOff = (shift > 0) ? shift : 0;

		const char* seq1 = seq_ + recOff;
		const char* seq2 = lzSeq_ + lzOff;
		uint32 minLen = MIN(seqLen_ - recOff, lzSeqLen_ - lzOff);

		int32 cc = insertCost;
		for (uint32 i = 0; i < minLen && cc < result_.cost.cost; ++i)
		{
			cc += (seq1[i] != seq2[i]) * classifierParams.mismatchCost;
		}

		if (cc < result_.cost.cost)
		{
			result_.cost.cost = cc;
			result_.cost.noMismatches = (cc - insertCost == 0);
			result_.shift = shift;

			return true;
		}

		return false;
	}


private:
	void PrepareLzBuffer();

	// TODO: optimize - these ones are the bottleneck
	// (1)
	MatchResult FindBestLzMatch(const FastqRecord& rec_,
								int32 recMinPos_,
								int32 maxEncodeThreshold_,
								int32 minCostThreshold_ = 0);

	const ReadsClassifierParams classifierParams;
	const MinimizerParameters minimParams;

	std::array<char, FastqRecord::MaxSeqLen> dummySequence;
        std::vector<std::array<char, FastqRecord::MaxSeqLen> > sequenceBuffer;
	std::deque<LzMatch*> lzBuffer;
};

#endif // READSCLASSIFIER_H
