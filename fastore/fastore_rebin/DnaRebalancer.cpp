/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "DnaRebalancer.h"

#include <algorithm>
#include <queue>


// helper functions
//
bool IsSizeGreaterThan(const MatchNode* node_, uint32 minSize_ = 0)
{
	uint64 size = 1;
	if (!node_->HasChildren())
		return size > minSize_;

	std::deque<MatchNode*> q(node_->children->begin(), node_->children->end());
	while (!q.empty())
	{
		auto n = q.front();
		q.pop_front();

		if (n->HasChildren())
		{
			size += n->children->size();
			q.insert(q.end(), n->children->begin(), n->children->end());
		}

		if (size > minSize_)
			break;
	}

	return size > minSize_;
}


void SwapNodeHierarchy(MatchNode* node_, MatchNode* prevNode_)
{
	if (node_->parentNode != NULL)
	{
		MatchNode* curParent = node_->parentNode;

		SwapNodeHierarchy(curParent, node_);

		curParent->RemoveChild(node_);
		curParent->parentNode = node_;		// INFO: lzRecord change too, but we don't need it atm
		node_->AddChild(curParent);
		node_->parentNode = prevNode_;
	}
	else
	{
		ASSERT(node_->type == MatchNode::TYPE_HARD);
		node_->type = MatchNode::TYPE_LZ;
	}
}


void SetAsRoot(MatchNode* node_)
{
	SwapNodeHierarchy(node_, NULL);
	node_->type = MatchNode::TYPE_HARD;
}


DnaRebalancer::DnaRebalancer(const MinimizerParameters& params_,
							 const BinBalanceParameters& binParams_,
							 bool pairedEnd_)
	:	FastqCategorizerBase(params_)
	,	binParams(binParams_)
	,	pairedEnd(pairedEnd_)
	,	readsClassifier(params_, binParams_.classifier)
{
	ASSERT(binParams_.validBinSignatures.size() > 1);
}


void DnaRebalancer::Rebalance(RebinContext& rebinCtx_,
								std::map<uint32, MatchNodesPtrBin>& bins_,
								uint32 signature_)
{
	ASSERT(rebinCtx_.graph->nodes.size() > 0);
	ASSERT(signature_ % binParams.signatureParity != 0);


	// TODO: optimization: reuse
	bins_.clear();

	// firstly, sort the reads
	//
	TFastqComparator<const MatchNode&> comparator;
	std::sort(rebinCtx_.graph->nodes.begin(), rebinCtx_.graph->nodes.end(), comparator);


	// build the tree to find HRs and LZs
	//
	readsClassifier.ConstructMatchTree(*rebinCtx_.graph, rebinCtx_.rootNodes);


	// traverse the reads graph
	//
	for (auto rootNode : rebinCtx_.rootNodes)
	{
		if (binParams.minTreeSize == 0 || (rootNode->HasChildren() && IsSizeGreaterThan(rootNode, binParams.minTreeSize)))
		{
			StoreTree(rootNode, rebinCtx_, bins_, signature_);
		}
		else
		{
			std::deque<MatchNode*> q(1, rootNode);
			do
			{
				MatchNode* n = q.front();
				q.pop_front();

				StoreSingleReadNode(n, bins_, signature_);

				if (n->HasChildren())
				{
					// WARN: after adding children to queue, we need to explicitely remove them
					// for compatibility with FastqNodePacker packing method,
					// otherwise we'll have stored multiple nodes duplicates
					q.insert(q.end(), n->children->begin(), n->children->end());
					n->RemoveChildren();
				}
			}
			while (!q.empty());
		}
	}
}



void DnaRebalancer::StoreTree(MatchNode* node_,
								RebinContext& rebinCtx_,
								std::map<uint32, MatchNodesPtrBin>& bins_,
								uint32 signature_)
{
	ASSERT(!node_->HasContigGroup());

	MatchNode* newRoot = node_;


	// find the new root(s)
	//
	std::tuple<uint32, uint16, bool> minimizerDesc = std::make_tuple(params.SignatureN(), 0, false);
	FastqRecordBuffer rcRec;
	rcRec.seqLen = node_->record->seqLen;

	if (binParams.selectMaxEdgeRead && newRoot->HasChildren())
	{
		// find the max left and min right node
		//
		const uint32 maxLeftPos = node_->record->seqLen - params.signatureLen - 1;
		const uint32 minRightPos = 0;

		auto leftRoot = std::make_pair(node_->record->minimPos, node_);
		auto rightRoot = leftRoot;

		std::deque<MatchNode*> nodesQueue(node_->children->begin(), node_->children->end());
		while (!nodesQueue.empty())
		{
			auto n = nodesQueue.front();
			nodesQueue.pop_front();

			if (leftRoot.first < n->record->minimPos)
				leftRoot = std::make_pair(n->record->minimPos, n);

			if (rightRoot.first > n->record->minimPos)
				rightRoot = std::make_pair(n->record->minimPos, n);

			if (leftRoot.first == maxLeftPos && rightRoot.first == minRightPos)
				break;

			if (n->HasChildren())
				nodesQueue.insert(nodesQueue.end(), n->children->begin(), n->children->end());
		}


		// now calculate the best signature for them
		//
		FastqRecordBuffer rcRecFwd, rcRecRev;

		std::tuple<uint32, uint16, bool> m1, m2;
		m1 = m2 = std::make_tuple(params.SignatureN(), 0, false);

		if (leftRoot.second != node_)
		{
			rcRecFwd.seqLen = node_->record->seqLen;
			leftRoot.second->record->ComputeRC(rcRecFwd);

			if (!pairedEnd)
				m1 = FindNewMinimizer(*leftRoot.second->record, rcRecFwd, signature_);
			else
			{
				auto mm = FindMinimizerHR(*leftRoot.second->record, signature_, binParams.signatureParity);
				m1 = std::make_tuple(mm.first, mm.second, false);
			}
		}

		if (rightRoot.second != node_)
		{
			rcRecRev.seqLen = node_->record->seqLen;
			rightRoot.second->record->ComputeRC(rcRecRev);

			if (!pairedEnd)
				m2 = FindNewMinimizer(*rightRoot.second->record, rcRecRev, signature_);
			else
			{
				auto mm = FindMinimizerHR(*rightRoot.second->record, signature_, binParams.signatureParity);
				m2 = std::make_tuple(mm.first, mm.second, false);
			}
		}


		// pick the best one
		//
		if (std::get<0>(m1) < std::get<0>(m2))
		{
			minimizerDesc = m1;
			newRoot = leftRoot.second;
			rcRec = std::move(rcRecFwd);
		}
		else if (std::get<0>(m2) != params.SignatureN())
		{
			minimizerDesc = m2;
			newRoot = rightRoot.second;
			rcRec = std::move(rcRecRev);
		}


		// swap the roots
		//
		if (newRoot != node_)
		{
			SetAsRoot(newRoot);
		}
	}
	else
	{
		// find the new signature for HR
		//
		FastqRecord* mainRec = node_->record;
		uint32 revSignaturePos = mainRec->seqLen - mainRec->minimPos - params.signatureLen;

		bool directionChange = false;

		rcRec.seqLen = mainRec->seqLen;
		mainRec->ComputeRC(rcRec);
		rcRec.minimPos = revSignaturePos;

		auto minimizerFwd = FindMinimizerHR(*mainRec, signature_, binParams.signatureParity);
		auto minimizerRev = !pairedEnd ? FindMinimizerHR(rcRec, signature_, binParams.signatureParity) : minimizerFwd;


		// temporarily do not search for rev-compl if the read has consensuses
		//
		std::pair<uint32, uint16> minimizer = minimizerFwd;
		if (minimizerFwd.first > minimizerRev.first)
		{
		   minimizer = minimizerRev;
		   directionChange = true;
		}

		minimizerDesc = std::make_tuple(minimizer.first, minimizer.second, directionChange);
	}

	FastqRecord* mainRec = newRoot->record;
	uint32 oldSignaturePos = mainRec->minimPos;


	// store record to bin
	//
	bool directionChange = std::get<2>(minimizerDesc);
	std::pair<uint32, uint16> minimizer = std::make_pair(std::get<0>(minimizerDesc), std::get<1>(minimizerDesc));

	// WARN: change of direction not supported in pared-end mode
	ASSERT(!directionChange || !pairedEnd);

	MatchNodesPtrBin* rb = NULL;
	if (minimizer.first != nBinValue)
	{
		ASSERT(minimizer.first != 0);
		rb = &bins_[minimizer.first];
		mainRec->minimPos = minimizer.second;
	}
	else	// we did not find any minimizer -- so use the current one
	{
		ASSERT(signature_ != 0);
		directionChange = false;		// no change has been made
		rb = &bins_[signature_];
	}

	// store the main encoding node
	//
	rb->nodes.push_back(newRoot);


	// WARN: do we need stats?
	//
	rb->stats.Update(*mainRec);


	// apply changes to the main record
	//
	if (directionChange)
	{
		mainRec->SetReadReverse(!mainRec->IsReadReverse());
		mainRec->CopyFrom(rcRec);

		// here we have already changed minimizer position
		oldSignaturePos = mainRec->seqLen - oldSignaturePos - params.signatureLen;
	}


	// check whether the direction change and new minimizer position affects record had exact matches
	//
	if (minimizer.first != nBinValue)
	{
		if (newRoot->HasExactMatches())
			UpdateExactMatches(newRoot, directionChange);

		if (newRoot->HasSubTreeGroup() && directionChange)
			UpdateTreeReads(newRoot, directionChange);
	}

	if (!newRoot->HasChildren())
		return;


	// update the children
	//
	uint64 treeSize = 0;
	std::deque<MatchNode*> q(newRoot->children->begin(), newRoot->children->end());
	while (!q.empty())
	{
		MatchNode* curNode = q.front();
		q.pop_front();
		treeSize++;

		// update the read itself
		//
		if (directionChange)
		{
			FastqRecord* curRec = curNode->record;
			curRec->ComputeRC(rcRec);
			curRec->SetReadReverse(!curRec->IsReadReverse());
			curRec->CopyFrom(rcRec);

			// update the signature position [!!!]
			//
			curRec->minimPos = curRec->seqLen - curRec->minimPos - params.signatureLen;
		}


		// update exact matches (if present)
		//
		if (curNode->HasExactMatches() && directionChange)
			UpdateExactMatches(curNode, directionChange);


		// update sub trees (if present)
		//
		if (curNode->HasSubTreeGroup() && directionChange)
			UpdateTreeReads(curNode, directionChange);


		// add children to processing queue
		//
		if (curNode->HasChildren())
		{
			q.insert(q.end(), curNode->children->begin(), curNode->children->end());
		}
	}

	// store the information about transferring of the tree
	//
	TreeTransferDefinition* g = rebinCtx_.CreateTransTreeGroup();
	g->signatureId = !directionChange ? signature_ : params.ReverseSignature(signature_);
	g->mainSignaturePos = oldSignaturePos;
	g->recordsCount = treeSize;
	newRoot->AddTransTreeGroup(g);
}


void DnaRebalancer::StoreSingleReadNode(MatchNode* node_, std::map<uint32, MatchNodesPtrBin>& bins_,
										  uint32 signature_)
{
	// find and select minimizers
	//
	// TODO: here we will operate directly on the indices
	//
	FastqRecord* rec = node_->record;
	const uint32 revSignaturePos = rec->seqLen - rec->minimPos - params.signatureLen;


	FastqRecordBuffer rcRec;
	rec->ComputeRC(rcRec);
	rcRec.minimPos = revSignaturePos;

	const bool allowRev = !pairedEnd || (!node_->HasExactMatches() && !node_->HasSubTreeGroup());

	auto minimizerFwd = FindMinimizerHR(*rec, signature_, binParams.signatureParity);
	auto minimizerRev = allowRev ? FindMinimizerHR(rcRec, signature_, binParams.signatureParity) : minimizerFwd;
	auto minimizer = minimizerFwd;

	bool directionChange = false;
	if (minimizerFwd.first > minimizerRev.first)
	{
	   minimizer = minimizerRev;
	   directionChange = true;
	}


	// store record to bin
	//
	MatchNodesPtrBin* rb = NULL;
	if (minimizer.first != nBinValue)
	{
		rb = &bins_[minimizer.first];
		rec->minimPos = minimizer.second;
	}
	else	// we did not find any minimizer -- so use the current one
	{
		directionChange = false;		// no change has been made
		rb = &bins_[signature_];
	}
	rb->nodes.push_back(node_);


	// WARN: do we need this???
	//
	rb->stats.Update(*rec);



	// apply changes to the main record
	//
	if (directionChange)
	{
		rec->SetReadReverse(!rec->IsReadReverse());
		rec->CopyFrom(rcRec);
	}


	// check whether the direction change and new minimizer position affects record had exact matches
	//
	if (minimizer.first != nBinValue)
	{
		if (node_->HasExactMatches())
			UpdateExactMatches(node_, directionChange);

		if (node_->HasSubTreeGroup() && directionChange)
			UpdateTreeReads(node_, directionChange);
	}
}


void DnaRebalancer::UpdateExactMatches(MatchNode* node_, bool directionChange_)
{
	const FastqRecord* mainRec = node_->record;
	ExactMatchesGroup* emg = node_->GetExactMatches();
	FastqRecordBuffer rcBuf;

	for (FastqRecord* em : emg->records)
	{
		// update the sequence and direction if changed
		//
		if (directionChange_)
		{
			// WARN: stanadrd copy does not handle PE nor quality
			//std::copy(mainRec->seq, mainRec->seq + mainRec->seqLen, em->seq);

			em->ComputeRC(rcBuf);
			em->CopyFrom(rcBuf);
			em->SetReadReverse(!em->IsReadReverse());
		}

		ASSERT(std::equal(mainRec->seq, mainRec->seq + mainRec->seqLen, em->seq));

		// update the signature position
		//
		em->minimPos = mainRec->minimPos;
	}
}


void DnaRebalancer::UpdateTreeReads(MatchNode* node_, bool directionChange_)
{
	ASSERT(node_->HasSubTreeGroup());

	if (!directionChange_)		// we need to change only the reads on reverse-compliment
		return;					// of the main encoding read


	const FastqRecord* mainRec = node_->record;

	FastqRecordBuffer rcRec;

	// when changind direction -- rev-compl, we also need to:
	// - rev-impl the sequence and signature
	// - inverse the ranges and variant positions
	//auto& groupList = reads_.trees.at(rh);

	//for (auto& tree : groupList)
	auto treeList = node_->GetSubTrees();
	for (auto tree: treeList)
	{
		// rev-compl the signature
		//
		uint32 sigRevLookup[4] = {3, 2, 1, 0};
		uint32 oldSig = tree->signatureId;
		tree->signatureId = 0;
		for (uint32 i = 0; i < params.signatureLen; ++i)
		{
			tree->signatureId <<= 2;
			tree->signatureId |= sigRevLookup[oldSig & 3];
			oldSig >>= 2;
		}

		// reverse tree signature position
		tree->mainSignaturePos = mainRec->seqLen - tree->mainSignaturePos - params.signatureLen;		// be careful here !!!


		std::array<char, MAX_SIGNATURE_LEN> newSignature;
		params.GenerateMinimizer(tree->signatureId, newSignature.data());

		ASSERT(std::search(mainRec->seq + tree->mainSignaturePos,
						   mainRec->seq + tree->mainSignaturePos + params.signatureLen,
						   newSignature.data(),
						   newSignature.data() + params.signatureLen) == mainRec->seq + tree->mainSignaturePos);


		// update the reads
		//
		for (MatchNode& mn : tree->nodes)
		{
			FastqRecord* rec = mn.record;

			// update the sequence
			//
			rec->ComputeRC(rcRec);
			rec->CopyFrom(rcRec);
			rec->SetReadReverse(!rec->IsReadReverse());

			// update the signature position
			//
			rec->minimPos = rec->seqLen - rec->minimPos - params.signatureLen;

			ASSERT(std::search(rec->seq + rec->minimPos,
							   rec->seq + rec->minimPos + params.signatureLen,
							   newSignature.data(),
							   newSignature.data() + params.signatureLen) == rec->seq + rec->minimPos);

			if (mn.HasExactMatches())
				UpdateExactMatches(&mn, true);

			if (mn.HasSubTreeGroup())
				UpdateTreeReads(&mn, true);
		}
	}
}


std::pair<uint32, uint16> DnaRebalancer::FindMinimizerHR(const FastqRecord &rec_,
														   uint32 curSignature,
														   uint32 curDivisor_)
{
	ASSERT(curDivisor_ > 1 && (curDivisor_ & (curDivisor_ - 1)) == 0);
	ASSERT(rec_.seqLen >= params.signatureLen - params.skipZoneLen);

	uint32 minimizer = maxLongMinimValue;
	uint32 pos = 0;

	for (int32 i = 0; i < rec_.seqLen - params.signatureLen - params.skipZoneLen; ++i)
	{
		const uint32 m = ComputeMinimizer(rec_.seq + i, params.signatureLen);

		if (m < minimizer &&
				m != curSignature &&
				m % curDivisor_ == 0 &&
				binParams.validBinSignatures[m] &&			// TODO: merge with valid signatures from categorizer
				IsMinimizerValid(m, params.signatureLen))
		{
			minimizer = m;
			pos = i;
		}
	}

	// filter the reads for which we cannot find a proper minimizer
	// and the ones which contain too much N symbols
	if (minimizer >= maxLongMinimValue
			|| std::count(rec_.seq, rec_.seq + rec_.seqLen, 'N') >= rec_.seqLen / 3)
		return std::make_pair(nBinValue, 0);

	return std::make_pair(minimizer, pos);
}


std::tuple<uint32, uint16, bool> DnaRebalancer::FindNewMinimizer(const FastqRecord &recFwd_,
																   const FastqRecord &recRev_,
																   uint32 curSignature_)
{
	auto minimizer = FindMinimizerHR(recFwd_, curSignature_, binParams.signatureParity);
	auto minimizerRev = FindMinimizerHR(recRev_,curSignature_, binParams.signatureParity);

	if (minimizer.first > minimizerRev.first)
		return std::make_tuple(minimizerRev.first, minimizerRev.second, true);

	return std::make_tuple(minimizer.first, minimizer.second, false);
}
