/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "ContigBuilder.h"
#include <queue>
#include <set>

#include <iostream>


ContigBuilder::ContigBuilder(const ReadsContigBuilderParams& params_,
								   const MinimizerParameters& minParams_)
	:	params(params_)
	,	minParams(minParams_)
{}

ContigBuilder::~ContigBuilder()
{}


void AddChildrenToQueue(std::deque<MatchNode*>& queue_, MatchNode* node_)
{
	ASSERT(node_->HasChildren());

	if (node_->children->size() > 1)
	{
		// prioritize
		std::vector<MatchNode*> cnodes;
		for (auto child : *node_->children)
		{
			if (!child->HasChildren())
				queue_.push_back(child);
			else
				cnodes.push_back(child);
		}
		queue_.insert(queue_.end(), cnodes.begin(), cnodes.end());
	}
	else
	{
		queue_.push_back(*node_->children->begin());
	}
}


bool ContigBuilder::Build(MatchNode* root_, PackContext& packCtx_)
{
	const uint32 seqLen = root_->record->seqLen;

	std::deque<MatchNode*> nextQueue;

	// do not encode the root node inside consesus
	//
	if (root_->HasChildren())
		AddChildrenToQueue(nextQueue, root_);


	// we shall keep the working nodes only here
	//
	ContigBuildInfo buildInfo;

	while (!nextQueue.empty())
	{
		// add the first node to consensus
		//
		MatchNode* node = nextQueue.front();
		nextQueue.pop_front();

		ASSERT(node->record->seqLen == seqLen);


		buildInfo.Reset(seqLen);


		// handle the case of invalid record not being added to the consensus
		// TODO: try to re-add the node later to the
		//
		if (!AddRecord(buildInfo, node))
		{
			if (node->HasChildren())
				AddChildrenToQueue(nextQueue, node);

			continue;
		}


		// proceed with building the consensus
		//
		std::deque<MatchNode*> curQueue = std::move(nextQueue);

		if (node->HasChildren())
			AddChildrenToQueue(curQueue, node);


		// firstly try to add nodes which fully match with conesus (no variants)
		//
		while (curQueue.size() > 0)
		{
			node = curQueue.front();
			curQueue.pop_front();

			if (!AddRecord(buildInfo, node, true))
			{
				nextQueue.push_back(node);
			}
			else
			{
				if (node->HasChildren())
					AddChildrenToQueue(curQueue, node);
			}
		}


		// now try to add nodes which can introduce variants, but are still valid
		//
		std::swap(curQueue, nextQueue);
		while (curQueue.size() > 0)
		{
			node = curQueue.front();
			curQueue.pop_front();

			if (!AddRecord(buildInfo, node))				// *
			{
				nextQueue.push_back(node);
			}
			else
			{
				if (node->HasChildren())
					AddChildrenToQueue(curQueue, node);
			}
		}


		// if consensus is invalid, then discard changes
		//
		if (buildInfo.nodes.size() < params.minConsensusSize)
			continue;


		// post-process the consensus
		//
		OptimizeContig(buildInfo);

		if (buildInfo.nodes.size() < params.minConsensusSize)
			continue;

		UpdateContigLinkage(buildInfo);

		PostProcessContig(buildInfo);

		StoreContig(buildInfo, packCtx_);
	}

	return true;
}


bool ContigBuilder::AddRecord(ContigBuildInfo& contig_, MatchNode* node_, bool tryFullMatchOnly_)
{
	ASSERT(node_->record->seqLen == contig_.consensus.readLen);
	ASSERT(node_->type == MatchNode::TYPE_LZ);

	if (AvoidTreesInConsensus && node_->HasSubTreeGroup())
		return false;

	const uint32 seqLen = contig_.consensus.readLen;
	const FastqRecord* rec = node_->record;
	const uint32 consBegin = seqLen - node_->record->minimPos;
	const uint32 consEnd = consBegin + seqLen;


	// try adding the node to already exising ones
	//
	if (contig_.nodes.size() > 0)
	{
		// calculate the Hamming distance and possibly induced new variants
		//
		ContigBuildInfo::WorkNode workNode;
		uint32 hammingDistance = 0;

		for (uint32 i = params.beginCut; i < seqLen - params.endCut; ++i)
		{
			uint32 p = consBegin + i;
			if (contig_.consensus.sequence[p] != ContigBuildInfo::EmptyChar && contig_.consensus.sequence[p] != rec->seq[i])
			{
				hammingDistance++;

				if (contig_.variantFreqPerPos[p] == 0)
					workNode.newVariantPositions.push_back(p);
			}
			// forbid adding Ns to consensus
			else if(contig_.consensus.sequence[p] == ContigBuildInfo::EmptyChar && rec->seq[i] == 'N')
			{
				return false;
			}
		}


		// check whether read qualifies to be added into consensus
		//
		if (tryFullMatchOnly_)
		{
			if (workNode.newVariantPositions.size() > 0)
				return false;
		}
		else
		{
			// TODO: improve the check
			//
			const uint32 maxShift = (params.maxRecordShiftDifference == 0) ? (seqLen / 2) : (params.maxRecordShiftDifference);

			if (!(workNode.newVariantPositions.size() == 0
				  || (hammingDistance <= params.maxHammingDistance
					  && workNode.newVariantPositions.size() <= params.maxNewVariantsPerRead))
				|| ABS((int32)contig_.nodes.back().match->record->minimPos - node_->record->minimPos) > (int32)maxShift
				|| GetConsEncodeCost(contig_, workNode, hammingDistance) > GetNormalEncodeCost(*node_))
			{
				return false;
			}
		}


		// introduce the new variants
		//
		contig_.consensus.variantsCount += workNode.newVariantPositions.size();


		// add the read to the consensus and update the stats
		//
		for (uint32 i = params.beginCut; i < seqLen - params.endCut; ++i)
		{
			uint32 p = consBegin + i;
			if (contig_.consensus.sequence[p] == ContigBuildInfo::EmptyChar)	// extend the consensus sequence
			{
				contig_.consensus.sequence[p] = rec->seq[i];
			}
			else if (contig_.consensus.sequence[p] != rec->seq[i])
			{
				contig_.variantFreqPerPos[p]++;
			}

			contig_.recordsPerPos[p]++;
		}

		// update the consensus range and add the record
		//

		// handle the case where the position of the signature is in the
		// initial padding range
		if (rec->minimPos <= params.beginCut)
			contig_.consensus.range.first = MIN(contig_.consensus.range.first, consBegin + rec->minimPos + minParams.signatureLen);
		else
			contig_.consensus.range.first = MIN(contig_.consensus.range.first, consBegin + params.beginCut);

		// handle the case where the position of the signature is in the
		// ending padding range
		if (seqLen - rec->minimPos - minParams.signatureLen <= params.endCut)
			contig_.consensus.range.second = MAX(contig_.consensus.range.second, consBegin + rec->minimPos);
		else
			contig_.consensus.range.second = MAX(contig_.consensus.range.second, consEnd - params.endCut);

		workNode.match = node_;
		contig_.nodes.push_back(std::move(workNode));
	}
	else
	{
		// discard reads with N's as the main consensus read						// TODO: remove this one
		//
		if (std::find(rec->seq, rec->seq + seqLen, 'N') != rec->seq + seqLen)		// OPT: use strchr()
			return false;

		// just insert the record into consensus as it is
		//
		for (uint32 i = params.beginCut; i < seqLen - params.endCut; ++i)			// OPT: use strcpy
			contig_.consensus.sequence[consBegin + i] = rec->seq[i];

		// update the consensus range and add the record
		//
		// handle the case where the position of the signature is in the
		// initial padding range
		if (rec->minimPos <= params.beginCut)
			contig_.consensus.range.first = consBegin + minParams.signatureLen + rec->minimPos;
		else
			contig_.consensus.range.first = consBegin + params.beginCut;

		// handle the case where the position of the signature is in the
		// ending padding range
		if (seqLen - rec->minimPos - minParams.signatureLen <= params.endCut)
			contig_.consensus.range.second = consBegin + rec->minimPos;
		else
			contig_.consensus.range.second = consEnd - params.endCut;

		ContigBuildInfo::WorkNode n;
		n.match = node_;
		contig_.nodes.push_back(std::move(n));
	}

	return true;
}


// TODO: remove removedNodes_ and use contig_ built-in
// TODO: refactor, remove 'idxToRemove'
uint32 ContigBuilder::OptimizeContig(ContigBuildInfo& contig_)
{
	// check whether we need to remove some records
	std::vector<uint32> idxToRemove;
	for (uint32 i = 0; i < contig_.nodes.size(); ++i)
	{
		ContigBuildInfo::WorkNode& n = contig_.nodes[i];
		if (n.newVariantPositions.size() > 0)
		{
			bool onlyOneVar = false;

			for (auto pos : n.newVariantPositions)
			{
				if (contig_.variantFreqPerPos[pos] == 1)
					onlyOneVar = true;
			}

			if (onlyOneVar)
				idxToRemove.push_back(i);			// TODO: remove the nodes here and use flag about removal
		}
	}


	if (idxToRemove.size() == 0)
		return 0;


	// filter the reads -- remove the records from the back
	auto oldNodes = std::move(contig_.nodes);

	contig_.Reset(contig_.consensus.readLen);
	uint32 idx = 0;
	uint32 removedNodesCount = 0;

	// TODO: refac
	for (std::vector<uint32>::iterator i = idxToRemove.begin(); i != idxToRemove.end(); ++i)
	{
		// CHECK
		if (*i > idx)
			contig_.nodes.insert(contig_.nodes.end(), oldNodes.begin() + idx, oldNodes.begin() + *i);

		contig_.removedNodes.push_back(oldNodes[*i].match);		// TODO: handle EMs
		removedNodesCount++;

		idx = *i + 1;
	}
	if (idx < oldNodes.size())
		contig_.nodes.insert(contig_.nodes.end(), oldNodes.begin() + idx, oldNodes.end());


	// rebuild the consensus
	//
	// set the first root record
	{
		const FastqRecord* rec = contig_.nodes[0].match->record;

		// just insert the record into consensus as it is
		const char* seq = rec->seq;
		uint32 consBegin = contig_.consensus.readLen - rec->minimPos;
		uint32 consEnd = consBegin + contig_.consensus.readLen;

		// use std::copy
		for (uint32 i = params.beginCut; i < contig_.consensus.readLen - params.endCut; ++i)
			contig_.consensus.sequence[consBegin + i] = seq[i];

		// handle the case where the position of the signature is in the
		// initial padding range
		if (rec->minimPos <= params.beginCut)
			contig_.consensus.range.first = consBegin + rec->minimPos + minParams.signatureLen;
		else
			contig_.consensus.range.first = consBegin + params.beginCut;

		// handle the case where the position of the signature is in the
		// ending padding range
		if (contig_.consensus.readLen - rec->minimPos - minParams.signatureLen <= params.endCut)
			contig_.consensus.range.second = consBegin + rec->minimPos;
		else
			contig_.consensus.range.second = consEnd - params.endCut;
	}


	// calculate the variants
	//
	for (std::vector<ContigBuildInfo::WorkNode>::iterator irn = contig_.nodes.begin() + 1; irn != contig_.nodes.end(); ++irn)
	{
		// calculate the possibly introduced new variants
		const FastqRecord* rec = irn->match->record;
		uint32 consBegin = contig_.consensus.readLen - rec->minimPos;
		uint32 consEnd = consBegin + contig_.consensus.readLen;

		irn->newVariantPositions.clear();

		for (uint32 i = params.beginCut; i < contig_.consensus.readLen - params.endCut; ++i)
		{
			uint32 p = consBegin + i;
			if (contig_.consensus.sequence[p] != ContigBuildInfo::EmptyChar && contig_.consensus.sequence[p] != rec->seq[i])
			{
				if (contig_.variantFreqPerPos[p] == 0)
					irn->newVariantPositions.push_back(p);
			}
		}

		// introduce the new variants
		contig_.consensus.variantsCount += irn->newVariantPositions.size();

		// add the read to the consensus and update the stats
		for (uint32 i = params.beginCut; i < contig_.consensus.readLen - params.endCut; ++i)
		{
			uint32 p = consBegin + i;
			if (contig_.consensus.sequence[p] == ContigBuildInfo::EmptyChar)	// extend the consensus sequence
				contig_.consensus.sequence[p] = rec->seq[i];
			else if (contig_.consensus.sequence[p] != rec->seq[i])
				contig_.variantFreqPerPos[p]++;

			contig_.recordsPerPos[p]++;
		}


		// update the consensus range and add the record
		//

		// handle the case where the position of the signature is in the
		// initial padding range
		if (rec->minimPos <= params.beginCut)
			contig_.consensus.range.first = MIN(contig_.consensus.range.first, consBegin + rec->minimPos + minParams.signatureLen);
		else
			contig_.consensus.range.first = MIN(contig_.consensus.range.first, consBegin + params.beginCut);

		// handle the case where the position of the signature is in the
		// ending padding range
		if (contig_.consensus.readLen - rec->minimPos - minParams.signatureLen <= params.endCut)
			contig_.consensus.range.second = MAX(contig_.consensus.range.second, consBegin + rec->minimPos);
		else
			contig_.consensus.range.second = MAX(contig_.consensus.range.second, consEnd - params.endCut);
	}

	return removedNodesCount;
}


void ContigBuilder::PostProcessContig(ContigBuildInfo &contig_)
{
	// average the conseusus
	//
	if (contig_.consensus.variantsCount != 0)
	{
		char dnaToIdx[128];
		std::fill(dnaToIdx, dnaToIdx + 128, -1);
		dnaToIdx[(int)'A'] = 0;
		dnaToIdx[(int)'G'] = 1;
		dnaToIdx[(int)'C'] = 2;
		dnaToIdx[(int)'T'] = 3;
		dnaToIdx[(int)'N'] = 4;
		const char* idxToDna = "AGCTN";

		// vote on average symbol
		//
		std::map<uint32, std::vector<uint32> > posStats;
		for (auto& workNode : contig_.nodes)
		{
			const char* seq = workNode.match->record->seq;
			uint32 consBegin = contig_.consensus.readLen - workNode.match->record->minimPos;

			for (uint32 i = params.beginCut; i < contig_.consensus.readLen - params.endCut; ++i)
			{
				uint32 p = consBegin + i;
				if (contig_.variantFreqPerPos[p] > 0)
				{
					if (posStats.count(p) == 0)
						posStats[p].resize(5);

					posStats[p].at(dnaToIdx[(int32)seq[i]]) += 1;
				}
			}
		}
		ASSERT(posStats.size() > 0);

		// update the consensus with the most frequent variants
		//
		for (auto& stat : posStats)
		{
			char c = 'N';
			uint32 maxFreq = 0 ;
			for (uint32 i = 0; i < 5; ++i)
			{
				if (stat.second[i] > maxFreq)
				{
					maxFreq = stat.second[i];
					c = idxToDna[i];
				}
			}
			//SOFT_ASSERT(c != 'N');					// we actually allow Ns
			contig_.consensus.sequence[stat.first] = c;
		}


		// update the variants count
		//
		for (uint32 i = 0; i < contig_.variantFreqPerPos.size(); ++i)
		{
			contig_.consensus.variantPositions[i] = contig_.variantFreqPerPos[i] != 0;
			contig_.consensus.variantsCount += contig_.variantFreqPerPos[i] != 0;
		}
	}

	// replace '.' symbols with 'N'
	//
	std::replace(contig_.consensus.sequence.begin(), contig_.consensus.sequence.end(),
				 (char)+ContigBuildInfo::EmptyChar, 'N');


	// sort by records minimizer position
	//
	std::sort(contig_.nodes.begin(), contig_.nodes.end());
}


void ContigBuilder::UpdateContigLinkage(ContigBuildInfo& contig_)
{
	// create node filters
	//

	std::vector<MatchNode*> consNodes;
	for (const auto& workNode : contig_.nodes)
		consNodes.push_back(workNode.match);
	std::sort(consNodes.begin(), consNodes.end());


	// search for parents outside the consensus and update them:
	// - remove the children
	// - select the only one parent (min LZ cost) and set consensus node as the only child
	//
	// WARNING: in the current situation we need to select the LAST parent, as all the parents
	// need to be encoded before encoding the consensus !!!
	//
	{
		MatchNode* bestParent = NULL;

		int32 minCost = 255;
		for (auto& n : contig_.nodes)
		{
			ASSERT(n.match->parentNode != NULL);

			if (!std::binary_search(consNodes.begin(), consNodes.end(), n.match->parentNode))
			{
				// unlink the contig nodes from parents which are outside contig
				//
				MatchNode* parent = n.match->parentNode;
				parent->RemoveChild(n.match);
				n.match->parentNode = NULL;

				contig_.unlinkedParents.push_back(parent);

				if (std::find(contig_.removedNodes.begin(), contig_.removedNodes.end(), parent) == contig_.removedNodes.end()
						&& (contig_.mainNode == NULL || minCost > n.match->encodeCost))
				{
					bestParent = parent;
					contig_.mainNode = n.match;
					minCost = n.match->encodeCost;
				}
			}
		}

		// link the consensus node with the rest of the tree (if not root)
		ASSERT(contig_.mainNode != NULL);
		ASSERT(bestParent != NULL);

		contig_.mainNode->parentNode = bestParent;
		bestParent->AddChild(contig_.mainNode);
	}

	// remove the main node from consensus nodes -- TODO: optimize
	//
	contig_.nodes.erase(std::find_if(contig_.nodes.begin(), contig_.nodes.end(),
									 [&](ContigBuildInfo::WorkNode& w_){return w_.match == contig_.mainNode;}));
	consNodes.erase(std::lower_bound(consNodes.begin(), consNodes.end(), contig_.mainNode));


	// remove children in main node if they reside in contig
	//
	if (contig_.mainNode->HasChildren())
	{
		for (auto ic = contig_.mainNode->children->begin(); ic != contig_.mainNode->children->end(); )
		{
			if (std::binary_search(consNodes.begin(), consNodes.end(), *ic))
				ic = contig_.mainNode->children->erase(ic);
			else
				ic++;
		}
	}


	// search for children outside the consensus and update them:
	// - set consensus node as parent
	// side effect: the lzRecord of the modified record differs from the record of the parent node
	//
	for (auto& n : contig_.nodes)
	{
		if (n.match->children == NULL)
			continue;

		for (auto ic = n.match->children->begin(); ic != n.match->children->end(); )
		{
			// is the child outside the concensus family?
			//
			if (!std::binary_search(consNodes.begin(), consNodes.end(), *ic))
			{
				// set the parent as consensus node
				(*ic)->parentNode = contig_.mainNode;
				contig_.mainNode->AddChild(*ic);

				ic = n.match->children->erase(ic);	// ***
			}
			else
			{
				ic++;		// ***: remember about earsing elements by iterators!!!
			}
		}


		// we can now also remove children, as they will point to
		// the nodes inside contig
		//
		n.match->RemoveChildren();
	}


	// handle the removed nodes:
	// - if parent not in consensus: leave untouched
	// - else: child of the consensus
	// side effect: the lzRecord of the modified record differs from the record of the parent node
	//
#if (DEV_DEBUG_MODE)
	for (MatchNode* n : contig_.removedNodes)
	{
		ASSERT(!std::binary_search(consNodes.begin(), consNodes.end(), n->parentNode));
	}
#endif
}


void ContigBuilder::StoreContig(ContigBuildInfo &contig_, PackContext &packCtx_)
{
	ASSERT(contig_.mainNode != NULL);
	ASSERT(contig_.nodes.size() > 0);

	ContigDefinition* contigDef = packCtx_.CreateContigGroup();
	contig_.mainNode->AddContigGroup(contigDef);

	// move the consensus info-  we won't need it later
	contigDef->consensus = std::move(contig_.consensus);


	// set the type of new reads and prepare the contig group
	//
	contig_.mainNode->type = MatchNode::TYPE_LZ;
	for (ContigBuildInfo::WorkNode& n : contig_.nodes)
	{
		n.match->type = MatchNode::TYPE_CONTIG_READ;
		contigDef->nodes.push_back(n.match);
	}
}
