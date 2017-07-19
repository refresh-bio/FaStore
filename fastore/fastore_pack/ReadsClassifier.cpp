/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "ReadsClassifier.h"

#include <algorithm>
#include <queue>
#include <set>
#include <array>
#include <functional>


ReadsClassifierSE::ReadsClassifierSE(const MinimizerParameters& minimParams_,
									 const ReadsClassifierParams& classifierParams_)
	:	classifierParams(classifierParams_)
	,	minimParams(minimParams_)
{
	std::fill(dummySequence.begin(), dummySequence.end(), 'N');

        sequenceBuffer.resize(classifierParams.maxLzWindowSize);
}


ReadsClassifierSE::~ReadsClassifierSE()
{
	for (uint32 i = 0; i < lzBuffer.size(); ++i)
		delete lzBuffer[i];
	lzBuffer.clear();
}


void ReadsClassifierSE::PrepareLzBuffer()
{
	for (uint32 i = 0; i < lzBuffer.size(); ++i)
		delete lzBuffer[i];
	lzBuffer.clear();

	for (uint32 i = 0; i < classifierParams.maxLzWindowSize; ++i)
	{
		LzMatch* lz = new LzMatch();
                std::copy(dummySequence.begin(), dummySequence.end(), sequenceBuffer[i].begin());
                lz->seqLen = sequenceBuffer[i].size();
				lz->seq = sequenceBuffer[i].data();
		lzBuffer.push_back(lz);
	}
}


ReadsClassifierSE::MatchResult ReadsClassifierSE::FindBestLzMatch(const FastqRecord &rec_,
																  int32 recMinPos_,
																  int32 maxCostThreshold_,
																  int32 minCostThreshold_)
{
	MatchResult result;
	result.cost.cost = maxCostThreshold_ + 1;

	for (uint32 i = 0; i < lzBuffer.size(); ++i)
	{
		const auto& lz = *lzBuffer[i];
		if (!UpdateLzMatchResult(result,
								 rec_.seq, rec_.seqLen, recMinPos_,
								 lz.seq, lz.seqLen, lz.minPos,
								 minCostThreshold_))
		{
			continue;
		}

		result.prevId = i;

		// found identical read or the cost is below threshold
		if (result.cost.cost == 0)
			break;
	}


	return result;
}

struct RootIdSorter
{
	bool operator() (const std::pair<const MatchNode*, uint32>& m1_,
					 const std::pair<const MatchNode*, uint32>& m2_) const
	{
		return m1_.second > m2_.second;
	}
};


void ReadsClassifierSE::ConstructMatchTree(GraphEncodingContext& graph_,
										   std::vector<MatchNode*>& outRootNodes_,
										   MatchNode* auxRootNode_)
{
	static const uint32 MinSignaturePos = 8;
	static const uint32 SigOffset = 2;
	static const uint32 BuffersPerPos = 5;

	outRootNodes_.clear();

	// take care of match nodes
	//
	PrepareLzBuffer();


	const bool usePrefixBuffer = classifierParams.extraReduceHardReads || classifierParams.extraReduceExpensiveLzMatches;


	// extra reverser-search buffer
	//
	auto prefixFun = [&](MatchNode* x_, MatchNode* y_)
	{

		ASSERT(x_->record->minimPos >= SigOffset);
		ASSERT(y_->record->minimPos >= SigOffset);

		const int32 maxRange = MIN(x_->record->minimPos, y_->record->minimPos) - SigOffset;

		const char* px = x_->record->seq + x_->record->minimPos - SigOffset;
		const char* py = y_->record->seq + y_->record->minimPos - SigOffset;

		for (int32 i = 0; i < maxRange; i++)
		{
			if (*px < *py)
				return true;
			if (*px > *py)
				return false;

			px--;
			py--;
		}

		return x_->record->minimPos > y_->record->minimPos;
	};

	typedef std::set<MatchNode*, std::function<bool(MatchNode*, MatchNode*)>> rp_buffer_t;

	std::array<rp_buffer_t*, BuffersPerPos*BuffersPerPos> rp_buffers;

	for (uint32 i = 0; i < rp_buffers.size(); ++i)
		rp_buffers[i] = usePrefixBuffer ? new rp_buffer_t(prefixFun) : NULL;

	int8 dnaToIdx[128];
	std::fill(dnaToIdx, dnaToIdx + 128, -1);
	dnaToIdx['A'] = 0;
	dnaToIdx['C'] = 1;
	dnaToIdx['G'] = 2;
	dnaToIdx['T'] = 3;
	dnaToIdx['N'] = 4;


	// handle the case when provided an auxilaty root
	//
	if (auxRootNode_ != NULL)
	{
		LzMatch* newLz = lzBuffer.back();
		lzBuffer.pop_back();

                // copy the sequence to local circular buffer to avoid L1 data cache misses
                // when searching for matches -- accounts for >70% miss when accessing the nested pointers
                //newLz->seq = auxRootNode_->record->seq;
                std::copy(auxRootNode_->record->seq, auxRootNode_->record->seq + auxRootNode_->record->seqLen, newLz->seq);

                newLz->seqLen = auxRootNode_->record->seqLen;
                newLz->node = auxRootNode_;
		newLz->minPos = auxRootNode_->record->minimPos;

		lzBuffer.push_front(newLz);

		outRootNodes_.push_back(auxRootNode_);
	}


	// perform matching
	//
	for (MatchNode& curNode : graph_.nodes)
	{
		FastqRecord* rec = curNode.record;

		LzMatch* newLz = lzBuffer.back();
		lzBuffer.pop_back();


		int32 encodeThreshold;
		if (classifierParams.encodeThresholdValue == 0)
			encodeThreshold = rec->seqLen / 2;				// automatic threshold
		else
			encodeThreshold = classifierParams.encodeThresholdValue;

		MatchResult matchResult = FindBestLzMatch(*rec, rec->minimPos, encodeThreshold);

		LzMatch* bestLz = lzBuffer[matchResult.prevId];

		// prepare new match
		//
				// copy the sequence to local circular buffer to avoid L1 data cache misses
				// when searching for matches -- accounts for >70% miss when accessing the nested pointers
				//newLz->seq = rec->seq;
				std::copy(rec->seq, rec->seq + rec->seqLen, newLz->seq);

				newLz->seqLen = rec->seqLen;
		newLz->node = &curNode;
		newLz->minPos = rec->minimPos;


		bool identicalReads = (matchResult.cost.cost == 0 && bestLz->seqLen == rec->seqLen);
		bool isHardRead = matchResult.cost.cost > encodeThreshold;


		// check whether parent is valid -- in case of compressing sub-trees
		//
		if (identicalReads)
		{
			identicalReads = lzBuffer[matchResult.prevId]->node->type != MatchNode::TYPE_NONE;
		}

		if (identicalReads)
		{
			curNode.type = MatchNode::TYPE_NONE;
			curNode.lzRecord = NULL;
			curNode.parentNode = NULL;


			// TODO: here we should only encode the fact that the node has exact matches
			// increasing the counter, not the number of nodes

			// WARN: not necesarily the last node will be the parent of the new LZ, as the
			// sorting starts from minimizer position and mismatches can be before that
			//uint32 lzNodeId = lzBuffer[matchResult.prevId]->nodeId;
			//MatchNode* parentNode = &graph_.nodes[lzNodeId];
			MatchNode* parentNode = lzBuffer[matchResult.prevId]->node;
			ASSERT(parentNode->type != MatchNode::TYPE_NONE);

			//FastqRecord* parentRec = parentNode->record;

#if DEV_DEBUG_MODE
			ASSERT(std::search(parentRec->seq,
							   parentRec->seq + parentRec->seqLen,
							   rec->seq, rec->seq + rec->seqLen)
					!= parentRec->seq + parentRec->seqLen);
#endif

			// move exact matches from this read to the main read
			//
			if (curNode.HasExactMatches())
			{
				ExactMatchesGroup* emg = curNode.GetExactMatches();

				if (!parentNode->HasExactMatches())
				{
					parentNode->CreateExactMatchesGroup(emg);
				}
				else
				{
					// merge exact matches groups
#if DEV_DEBUG_MODE
					for (FastqRecord* r : emg->records)
						ASSERT(r->minimPos == parentRec->minimPos);
#endif

					ExactMatchesGroup* parentEmg = parentNode->GetExactMatches();
					parentEmg->records.insert(parentEmg->records.end(), emg->records.begin(), emg->records.end());
				}

				curNode.RemoveExactMatches();
			}

			// move trees from the current read to the parentread
			//
			if (curNode.HasSubTreeGroup())
			{
				auto subTrees = curNode.GetSubTrees();

				for (auto tree : subTrees)
				{
					// here we can optimize to operate directly on decorators
					parentNode->AddSubTreeGroup((GraphEncodingContext*)tree);
				}

				curNode.RemoveSubTrees();
			}


			// add the exact match itself to parent's node
			//
			ExactMatchesGroup* parentEmg = NULL;
			if (!parentNode->HasExactMatches())
			{
				parentEmg = graph_.CreateExactMatchesGroup();
				parentNode->CreateExactMatchesGroup(parentEmg);
			}
			else
			{
				parentEmg = parentNode->GetExactMatches();
			}
			parentEmg->records.push_back(rec);


			// update lz buffer
			//
			lzBuffer.push_back(newLz);
		}
		else
		{
			MatchNode* parentNode = NULL;
			rp_buffer_t* rp_buffer = NULL;

			if (usePrefixBuffer)
			{
				// TODO: another option to tweak
				//
				const uint32 expensiveLzThreshold = encodeThreshold / 2;

				const bool searchInRevBuffer = isHardRead
						|| (classifierParams.extraReduceExpensiveLzMatches && matchResult.cost.cost > (int32)expensiveLzThreshold);

				// if we also search for better match for LZ
				//
				if (!isHardRead)
					encodeThreshold = expensiveLzThreshold;



				// assign the rp_buffer - we'll use it later to filter the appropriate reads to be added into
				if (curNode.record->minimPos >= MinSignaturePos)
				{
					int32 bufIdx = dnaToIdx[(int32)curNode.record->seq[curNode.record->minimPos - 2]];
					bufIdx = bufIdx * BuffersPerPos + dnaToIdx[(int32)curNode.record->seq[curNode.record->minimPos - 1]];

					ASSERT(bufIdx >= 0 && bufIdx < (int32)(BuffersPerPos * BuffersPerPos));
					rp_buffer = rp_buffers[bufIdx];
				}

				if (searchInRevBuffer && curNode.record->minimPos >= MinSignaturePos)
				{
					std::pair<MatchResult, MatchNode*> resultFwd, resultRev;
					resultFwd.first.cost.cost = encodeThreshold + 1;
					resultRev.first.cost.cost = encodeThreshold + 1;

					auto p = rp_buffer->lower_bound(&curNode);
					std::set<MatchNode*, std::function<bool(MatchNode*, MatchNode*)>>::reverse_iterator q(p);

					for (uint32 cnt = 0; p != rp_buffer->end() && cnt < classifierParams.maxLzWindowSize / 2 + 1; ++p, ++cnt)
					{
						// provisional LzMatch -- TODO: optimize
						const auto lz = (*p)->record;
						if (!UpdateLzMatchResult(resultFwd.first,
												 curNode.record->seq, curNode.record->seqLen, curNode.record->minimPos,
												 lz->seq, lz->seqLen, lz->minimPos))
						{
							continue;
						}

						resultFwd.second = *p;
					}

					for (uint32 cnt = 0; q != rp_buffer->rend() && cnt < classifierParams.maxLzWindowSize / 2 + 1; ++q, ++cnt)
					{
						const auto lz = (*q)->record;
						if (!UpdateLzMatchResult(resultRev.first,
												 curNode.record->seq, curNode.record->seqLen, curNode.record->minimPos,
												 lz->seq, lz->seqLen, lz->minimPos))
						{
							continue;
						}

						resultRev.second = *q;
					}

					const int32 minCost = MIN(resultFwd.first.cost.cost, resultRev.first.cost.cost);


					if (minCost < encodeThreshold && minCost < matchResult.cost.cost)
					{
						if (resultFwd.first.cost.cost < resultRev.first.cost.cost)
						{
							parentNode = resultFwd.second;
							matchResult = resultFwd.first;
						}
						else
						{
							parentNode = resultRev.second;
							matchResult = resultRev.first;
						}

						isHardRead = false;
					}
				}
			}


			// standard LZ?
			if (isHardRead)
			{
				// the read has no parent -- a possible root match, do not link with graph
				curNode.type = MatchNode::TYPE_HARD;
				curNode.lzRecord = NULL;
				curNode.parentNode = NULL;

				outRootNodes_.push_back(&curNode);
			}
			else
			{
				if (parentNode == NULL)
				{
					parentNode = lzBuffer[matchResult.prevId]->node;
				}


				// special case when linking to exact matches???
				ASSERT(matchResult.shift == parentNode->record->minimPos - curNode.record->minimPos);

				curNode.type = MatchNode::TYPE_LZ;
				curNode.parentNode = parentNode;
				curNode.lzRecord = parentNode->record;
				curNode.shiftValue = matchResult.shift;
				curNode.SetNoMismatches(matchResult.cost.noMismatches);
				curNode.encodeCost = matchResult.cost.cost;

				parentNode->AddChild(&curNode);
			}

			lzBuffer.push_front(newLz);

			if (rp_buffer != NULL)
				rp_buffer->insert(&curNode);
		}
	}


	// cleanup
	//
	for (auto rp : rp_buffers)
	{
		if (rp != NULL)
			delete rp;
	}
}

