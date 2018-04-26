/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"

#include <deque>

#include "NodesPacker.h"
#include "BitMemory.h"

#include "../fastore_bin/Params.h" // for BinModuleConfig definition


std::pair<uint32, uint32> GetSizeClassAndBits(uint64 groupSize_)
{
	ASSERT(groupSize_ > 0);

	uint32 emSizeClass = 0;
	uint32 bitsPerSize = 4;
	while (groupSize_ >= (1ULL << bitsPerSize))
	{
		emSizeClass++;
		bitsPerSize *= 2;
	}
	ASSERT(bitsPerSize <= 32);

	bitsPerSize = MIN(bitsPerSize, 30);

	ASSERT(groupSize_ <= (1ULL << bitsPerSize));

	return std::make_pair(emSizeClass, bitsPerSize);
}

uint32 GetBitsPerClass(uint32 classId_)
{
	ASSERT(classId_ < 4);

	const uint32 bitsPerClass[] = {4, 8, 16, 30};
	return bitsPerClass[classId_];
}


void IFastqNodesPacker::PackToBins(const std::map<uint32, MatchNodesPtrBin>& dnaBins_,
								  BinaryBinBlock& binBlock_)
{
	binBlock_.Clear();
	binBlock_.blockType = BinaryBinBlock::MultiSignatureType;

	BitMemoryWriter metaWriter(binBlock_.metaData);
	BitMemoryWriter dnaWriter(binBlock_.dnaData);
	BitMemoryWriter quaWriter(binBlock_.quaData);
	BitMemoryWriter headWriter(binBlock_.headData);


	// store standard bins
	//
	const uint32 nBinId = binConfig.minimizer.TotalMinimizersCount();
	const auto nBinIter = dnaBins_.find(nBinId);

	uint64 totalRecordsCount = 0;
	for (auto i = dnaBins_.cbegin(); i != nBinIter; ++i)
	{
		uint32 binId = i->first;
		const auto& bin = i->second;

		// skip empty bins
		ASSERT(bin.nodes.size() != 0);
		ASSERT(binId != 0);		// DEBUG

		BinaryBinDescriptor& desc = binBlock_.descriptors[binId];
		desc.Clear();

#if DEV_DEBUG_MODE
		char minString[64];
		binConfig.minimizer.GenerateMinimizer(binId, minString);

		for (const auto n : i->second.nodes)
		{
			const char* minSeq = std::search(n->record->seq,
											 n->record->seq + n->record->seqLen,
											 minString,
											 minString + binConfig.minimizer.signatureLen);

			ASSERT(minSeq != n->record->seq + n->record->seqLen);
			ASSERT(minSeq == n->record->seq + n->record->minimPos);
		}
#endif

		PackToBin(bin, metaWriter, dnaWriter, quaWriter, headWriter, desc, false);

		binBlock_.rawDnaSize += desc.rawDnaSize;
		binBlock_.rawHeadSize += desc.rawHeadSize;

		totalRecordsCount += desc.recordsCount;
	}


	// pack last N bin -- if present
	//
	if (nBinIter != dnaBins_.cend() && dnaBins_.at(nBinId).nodes.size() > 0)
	{
		BinaryBinDescriptor& nDesc = binBlock_.descriptors[nBinId];
		nDesc.Clear();

		const auto& bin = dnaBins_.at(nBinId);

		PackToBin(bin, metaWriter, dnaWriter, quaWriter, headWriter, nDesc, true);

		binBlock_.rawDnaSize += (uint64)nDesc.rawDnaSize;
		binBlock_.rawHeadSize += (uint64)nDesc.rawHeadSize;

		totalRecordsCount += (uint64)nDesc.recordsCount;
	}

	binBlock_.metaSize = metaWriter.Position();
	binBlock_.dnaSize = dnaWriter.Position();
	binBlock_.quaSize = quaWriter.Position();
	binBlock_.headSize = headWriter.Position();

	//ASSERT(totalRecordsCount == reads_.rawRecords.size());
}


void IFastqNodesPacker::PackToBin(const GraphEncodingContext& graph_,
								 BinaryBinBlock& binBin_,
								 uint32 signatureId_)
{
	binBin_.Clear();
	binBin_.blockType = BinaryBinBlock::SingleSignatureType;
	binBin_.signature = signatureId_;

	BitMemoryWriter metaWriter(binBin_.metaData);
	BitMemoryWriter dnaWriter(binBin_.dnaData);
	BitMemoryWriter quaWriter(binBin_.quaData);
	BitMemoryWriter headWriter(binBin_.headData);

	binBin_.auxDescriptors.resize(1);
	BinaryBinDescriptor& desc = binBin_.auxDescriptors.back();

	bool nBin = (signatureId_ == binConfig.minimizer.TotalMinimizersCount());

	// TODO: fix it!
	// it's a temporary temporary wrapper
	MatchNodesPtrBin matchBin;
	for (auto& node : graph_.nodes)
	{
		matchBin.nodes.push_back((MatchNode*)&node);
		matchBin.stats.Update(*node.record);
	}

	PackToBin(matchBin, metaWriter, dnaWriter, quaWriter, headWriter, desc, nBin);

	binBin_.rawDnaSize += (uint64)desc.rawDnaSize;
	binBin_.rawHeadSize += (uint64)desc.rawHeadSize;

	binBin_.metaSize = metaWriter.Position();
	binBin_.dnaSize = dnaWriter.Position();
	binBin_.quaSize = quaWriter.Position();
	binBin_.headSize = headWriter.Position();
}


// TODO: split into sub functions and templatize
//
void IFastqNodesPacker::PackToBin(const MatchNodesPtrBin& fqBin_,
								 BitMemoryWriter& metaWriter_,
								 BitMemoryWriter& dnaWriter_,
								 BitMemoryWriter& quaWriter_,
								 BitMemoryWriter& headWriter_,
								 BinaryBinDescriptor& binDesc_,
								 bool nBin_)
{
	const uint64 initialMetaPos = metaWriter_.Position();
	const uint64 initialDnaPos = dnaWriter_.Position();
	const uint64 initialQuaPos = quaWriter_.Position();
	const uint64 initialHeadPos = headWriter_.Position();


	binDesc_.recordsCount = 0;

	// verify the bin
	//
	ASSERT(fqBin_.nodes.size() != 0);
	ASSERT(fqBin_.nodes.size() < (1 << 28));
	ASSERT(fqBin_.stats.minSeqLen > 0);
	ASSERT(fqBin_.stats.maxSeqLen > 0);
	ASSERT(fqBin_.stats.maxSeqLen >= fqBin_.stats.minSeqLen);


	// initialize packing settings
	//
	BinPackSettings settings;
	settings.minLen = fqBin_.stats.minSeqLen;
	settings.maxLen = fqBin_.stats.maxSeqLen;
	settings.hasConstLen = (settings.minLen == settings.maxLen);
	settings.hasReadGroups = true;
	settings.usesHeaders = binConfig.archiveType.readsHaveHeaders;

	if (!nBin_)
		settings.suffixLen = binConfig.minimizer.signatureLen;
	else
		settings.suffixLen = 0;

	// if needed, the lenghts can be Huffman'd
	if (!settings.hasConstLen)
		settings.bitsPerLen = bit_length(fqBin_.stats.maxSeqLen - fqBin_.stats.minSeqLen);


	// store meta data header
	//
	metaWriter_.PutBits(fqBin_.stats.minSeqLen, LenBits);
	metaWriter_.PutBits(fqBin_.stats.maxSeqLen, LenBits);
	metaWriter_.PutBit(settings.hasReadGroups);


	// start packing
	//
	std::deque<const MatchNode*> packContext(fqBin_.nodes.begin(), fqBin_.nodes.end());

	while (!packContext.empty())
	{
		auto n = packContext.front();
		packContext.pop_front();

		StoreNextNode(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings, binDesc_, n, packContext);
	}


	// end packing
	//
	metaWriter_.FlushPartialWordBuffer();
	dnaWriter_.FlushPartialWordBuffer();
	quaWriter_.FlushPartialWordBuffer();
	headWriter_.FlushPartialWordBuffer();

	//ASSERT(metaWriter_.Position() - initialMetaPos < (1ULL << 32));
	//ASSERT(dnaWriter_.Position() - initialDnaPos < (1ULL << 32));
	//ASSERT(quaWriter_.Position() - initialQuaPos < (1ULL << 32));
	//ASSERT(headWriter_.Position() - initialHeadPos < (1ULL << 32));

	binDesc_.metaSize = metaWriter_.Position() - initialMetaPos;
	binDesc_.dnaSize = dnaWriter_.Position() - initialDnaPos;
	binDesc_.quaSize = quaWriter_.Position() - initialQuaPos;
	binDesc_.headSize = headWriter_.Position() - initialHeadPos;
}


void IFastqNodesPacker::StoreNextNode(BitMemoryWriter& metaWriter_,
									   BitMemoryWriter& dnaWriter_,
									  BitMemoryWriter& quaWriter_,
									  BitMemoryWriter& headWriter_,
									  const BinPackSettings& settings_,
									   BinaryBinDescriptor& binDesc_,
									   const MatchNode* node_,
									   std::deque<const MatchNode*>& packContext_)
{

	// store single read + extra metainfo
	//
	StoreRecordData(metaWriter_, dnaWriter_, quaWriter_, headWriter_,
					settings_, binDesc_, *node_->record);

	//ASSERT((uint64)binDesc_.recordsCount + 1 < (1UL << 32));
	binDesc_.recordsCount++;


	// store extra matching information
	//

	// store exact matches of the read (if present)
	//
	ASSERT(settings_.hasReadGroups);

	metaWriter_.PutBit(node_->HasExactMatches());
	metaWriter_.PutBit(node_->HasSubTreeGroup() || node_->HasTransTreeGroup());

	if (node_->HasExactMatches())
	{
		//const uint32 EmGroupSizeBits = 16;
		const ExactMatchesGroup* group = node_->GetExactMatches();
		ASSERT(group->records.size() > 0);

		auto classBits = GetSizeClassAndBits(group->records.size());
		ASSERT(group->records.size() < (1ULL << classBits.second));

		metaWriter_.Put2Bits(classBits.first);
		metaWriter_.PutBits(group->records.size(), classBits.second);


		for (const FastqRecord* em: group->records)
		{
			ASSERT(std::equal(node_->record->seq,
							  node_->record->seq + node_->record->seqLen,
							  em->seq));

			ASSERT(em->minimPos == node_->record->minimPos);

			StoreExactMatch(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, binDesc_, *em);

			//ASSERT((uint64)binDesc_.recordsCount + 1 < (1UL << 32));
			binDesc_.recordsCount++;
		}
	}


	//bool encodeChildren = true;
	if (node_->HasSubTreeGroup() || node_->HasTransTreeGroup())
	{
		std::vector<GraphEncodingContext*> treeList;

		// calculate and store the number of trees
		//
		uint64 treeCount = 0;

		if (node_->HasTransTreeGroup())
			treeCount++;

		if (node_->HasSubTreeGroup())
		{
			treeList = node_->GetSubTrees();
			treeCount += treeList.size();
		}

		auto classBits = GetSizeClassAndBits(treeCount);
		ASSERT(treeCount < (1ULL << classBits.second));

		metaWriter_.Put2Bits(classBits.first);
		metaWriter_.PutBits(treeCount, classBits.second);


		// store firstly the transfer tree
		//
		if (node_->HasTransTreeGroup())
		{
			auto tree = node_->GetTransTree();

			std::deque<const MatchNode*> tg(node_->children->begin(), node_->children->end());
			StoreNextGroup(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, binDesc_,
						   tree->signatureId, tree->mainSignaturePos, tree->recordsCount,
						   tg);
		}


		// store the sub trees
		//
		for (const GraphEncodingContext* group : treeList)
		{
			std::deque<const MatchNode*> tg;
			for (const auto& n : group->nodes)
			{
				tg.push_back(&n);
			}

			StoreNextGroup(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, binDesc_,
						   group->signatureId, group->mainSignaturePos, group->nodes.size(),
						   tg);
		}
	}


	// add the children
	//
	if (node_->HasChildren() && !node_->HasTransTreeGroup())
	{
		packContext_.insert(packContext_.end(), node_->children->begin(), node_->children->end());
	}
}


void IFastqNodesPacker::StoreNextGroup(BitMemoryWriter& metaWriter_,
									  BitMemoryWriter& dnaWriter_,
									   BitMemoryWriter& quaWriter_,
									   BitMemoryWriter& headWriter_,
									   const BinPackSettings& settings_,
									  BinaryBinDescriptor& binDesc_,
									  uint32 signatureId_,
									  uint32 signaturePos_,
									  uint64 groupSize_,
									  std::deque<const MatchNode*>& nodes_)
{
	ASSERT(nodes_.size() > 0);
	ASSERT(nodes_.size() <= groupSize_);

	// store the tree definition
	//
	const uint32 bitsPerSignature = binConfig.minimizer.signatureLen * 2;
	metaWriter_.PutBits(signatureId_, bitsPerSignature);
	metaWriter_.PutBits(signaturePos_, 8);

	auto classBits = GetSizeClassAndBits(groupSize_);
	ASSERT(groupSize_ < (1ULL << classBits.second));

	metaWriter_.Put2Bits(classBits.first);		// total records count in cur level
	metaWriter_.PutBits(groupSize_, classBits.second);		// total records count in cur level


#if(DEV_DEBUG_MODE)
	char sig[MAX_SIGNATURE_LEN];
	binConfig.minimizer.GenerateMinimizer(signatureId_, sig);
#endif

	// WARN: reusing the passed nodes deque
	//
	while (!nodes_.empty())
	{
		auto n = nodes_.front();
		nodes_.pop_front();

#if(DEV_DEBUG_MODE)
		ASSERT(std::search(n->record->seq, n->record->seq + n->record->seqLen,
						   sig, sig + binConfig.minimizer.signatureLen) != n->record->seq + n->record->seqLen);
#endif

		StoreNextNode(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, binDesc_, n, nodes_);
	}
}


void IFastqNodesPacker::UnpackFromBin(const BinaryBinBlock& binBin_,
									 std::vector<FastqRecord>& reads_,
									 GraphEncodingContext& graph_,
									 FastqRecordBinStats& stats_,
									 FastqChunk& fqChunk_,
									 bool append_)
{
	ASSERT(binBin_.blockType == BinaryBinBlock::SingleSignatureType);
	ASSERT(binBin_.auxDescriptors.size() > 0);
	ASSERT(binBin_.signature != 0);


	// calculate global bin properties
	//
	uint64 recordsCount = 0;

	for (const auto& desc : binBin_.auxDescriptors)
	{
		ASSERT(desc.recordsCount > 0);
		recordsCount += (uint64)desc.recordsCount;
	}

	ASSERT(recordsCount < (1 << 28));
	ASSERT(recordsCount != 0);

	// initialize buffers
	//
	uint64 recId = 0;
	if (append_)
	{
		recId = reads_.size();
		recordsCount += recId;

		const uint64 maxSize = binBin_.rawDnaSize * 2	// * 2 <-- to handle quality
				+ ((uint64)binConfig.archiveType.readsHaveHeaders * binBin_.rawHeadSize);

		ASSERT(fqChunk_.data.Size() >= fqChunk_.size + maxSize);
	}
	else
	{
		const uint64 chunkSize = binBin_.rawDnaSize * 2	// * 2 <-- to handle quality
				+ ((uint64)binConfig.archiveType.readsHaveHeaders * binBin_.rawHeadSize);
		if (fqChunk_.data.Size() < chunkSize)
			fqChunk_.data.Extend(chunkSize);
	}

	// reset stats for initialization
	//
	if (recId == 0)
	{
		ASSERT(fqChunk_.size == 0);
		ASSERT(reads_.size() == 0);

		stats_.Clear();
	}

	reads_.resize(recordsCount);


	// initialize IO
	//
	BitMemoryReader metaReader(binBin_.metaData, binBin_.metaSize);
	BitMemoryReader dnaReader(binBin_.dnaData, binBin_.dnaSize);
	BitMemoryReader quaReader(binBin_.quaData, binBin_.quaSize);
	BitMemoryReader headReader(binBin_.headData, binBin_.headSize);

	// initialize settings
	//
	BinPackSettings settings;
	settings.signatureId = binBin_.signature;

	if (binBin_.signature != binConfig.minimizer.TotalMinimizersCount())	// N bin?
	{
		settings.suffixLen = binConfig.minimizer.signatureLen;
		binConfig.minimizer.GenerateMinimizer(binBin_.signature, settings.signatureString.data());
	}
	else
	{
		settings.suffixLen = 0;
	}
	settings.usesHeaders = binConfig.archiveType.readsHaveHeaders;


	// unpack records
	//
	for (const BinaryBinDescriptor& desc : binBin_.auxDescriptors)
	{
		ASSERT(desc.recordsCount > 0);

		const uint64 initialMetaPos = metaReader.Position();
		const uint64 initialDnaPos = dnaReader.Position();
		const uint64 initialQuaPos = quaReader.Position();
		const uint64 initialHeadPos = headReader.Position();

		// read bin header
		//
		settings.minLen = metaReader.GetBits(LenBits);
		settings.maxLen = metaReader.GetBits(LenBits);
		ASSERT(settings.minLen > 0);
		ASSERT(settings.maxLen >= settings.minLen);

		settings.hasReadGroups = metaReader.GetBit() != 0;
		settings.hasConstLen = (settings.minLen == settings.maxLen);

		if (!settings.hasConstLen)
		{
			settings.bitsPerLen = bit_length(settings.maxLen - settings.minLen);
			ASSERT(settings.bitsPerLen > 0);
		}

		// read records from bin
		//
		uint64 totalReads = recId + (uint64)desc.recordsCount;

		const uint64 initialChunkSize = fqChunk_.size;
		while (recId < totalReads)
		{
			// todo: store the number of total nodes to avoid copying with resizing
			graph_.nodes.push_back(MatchNode());
			MatchNode& node = graph_.nodes.back();

			ReadNextNode(graph_, node,  reads_,  recId,
						 metaReader, dnaReader, quaReader, headReader,
						 settings, fqChunk_);
		}

		metaReader.FlushInputWordBuffer();
		dnaReader.FlushInputWordBuffer();
		quaReader.FlushInputWordBuffer();
		headReader.FlushInputWordBuffer();

		ASSERT(metaReader.Position() - initialMetaPos == (uint64)desc.metaSize);
		ASSERT(dnaReader.Position() - initialDnaPos == (uint64)desc.dnaSize);
		ASSERT(quaReader.Position() - initialQuaPos == (uint64)desc.quaSize);
		ASSERT(headReader.Position() - initialHeadPos == (uint64)desc.headSize);

		ASSERT(fqChunk_.size == initialChunkSize + (uint64)desc.rawDnaSize * 2 + (uint64)desc.rawHeadSize);
	}


	// update the stats
	//
	stats_.minSeqLen = settings.minLen;
	stats_.maxSeqLen = settings.maxLen;
}


void IFastqNodesPacker::ReadNextNode(GraphEncodingContext& graph_,
									MatchNode& curNode_,
									std::vector<FastqRecord>& reads_,
									uint64& recIdx_,
									BitMemoryReader& metaReader_,
									BitMemoryReader& dnaReader_,
									BitMemoryReader& quaReader_,
									BitMemoryReader& headReader_,
									const BinPackSettings& settings_,
									FastqChunk& fqChunk_)
{
	ASSERT(reads_.size() > recIdx_);

	// read the record data
	//
	FastqRecord& mainRec = reads_[recIdx_++];
	mainRec.Reset();

	ReadRecordData(metaReader_, dnaReader_, quaReader_, headReader_,
				   settings_, fqChunk_, mainRec);


	// bind the record to the current node
	//
	curNode_.record = &mainRec;



	// read the additional record groups (if present)
	//
	if (settings_.hasReadGroups)
	{
		// now it depends whether we should extract to read groups and consensuses
		// or to read array directly
		bool hasExactMatchesGroups = metaReader_.GetBit() != 0;
		bool hasTreeGroups = metaReader_.GetBit() != 0;

		if (hasExactMatchesGroups)
		{
			const uint32 emSizeClass = metaReader_.Get2Bits();
			const uint32 bits = GetBitsPerClass(emSizeClass);
			const uint32 groupSize = metaReader_.GetBits(bits);

			ASSERT(groupSize != 0);

			ExactMatchesGroup* emg = graph_.CreateExactMatchesGroup();
			//emg.records.resize(groupSize);
			curNode_.CreateExactMatchesGroup(emg);

			for (uint32 i = 0; i < groupSize; ++i)
			{
				ASSERT(settings_.maxLen == settings_.minLen);

				FastqRecord& em = reads_[recIdx_++];

				ReadExactMatch(metaReader_, dnaReader_, quaReader_, headReader_,
							   settings_, fqChunk_, mainRec, em);

				emg->records.push_back(&em);
			}
		}

		if (hasTreeGroups)
		{
			const uint32 tSizeClass = metaReader_.Get2Bits();
			const uint32 bits = GetBitsPerClass(tSizeClass);
			const uint32 tCount = metaReader_.GetBits(bits);

			ASSERT(tCount > 0);

			for (uint32 t = 0; t < tCount; ++t)
			{
				// read group meta information
				//
				GraphEncodingContext* group = graph_.CreateSubTreeGroup();
				curNode_.AddSubTreeGroup(group);

				// read consensus definition
				//
				const uint32 bitsPerSignature = binConfig.minimizer.signatureLen * 2;
				group->signatureId = metaReader_.GetBits(bitsPerSignature);
				group->mainSignaturePos = metaReader_.GetBits(8);

				const uint32 rcSizeClass = metaReader_.Get2Bits();
				const uint32 bits = GetBitsPerClass(rcSizeClass);
				const uint32 groupSize = metaReader_.GetBits(bits);

				ASSERT(groupSize != 0);

				// the reads inside the tree share different main signature
				// than the one from the beginning
				BinPackSettings consSettings = settings_;
				consSettings.signatureId = group->signatureId;
				consSettings.suffixLen = binConfig.minimizer.signatureLen;
				binConfig.minimizer.GenerateMinimizer(consSettings.signatureId, consSettings.signatureString.data());

				// read group
				//
				group->nodes.resize(groupSize);
				for (uint32 i = 0; i < groupSize; ++i)
				{
					MatchNode& n = group->nodes[i];

					// already read the contents of the subgroup
					//
					ReadNextNode(*group, n, reads_, recIdx_,
								 metaReader_, dnaReader_, quaReader_, headReader_,
								 consSettings, fqChunk_);
				}
			}
		}
	}
}


// SE packer
//
void FastqNodesPackerSE::StoreRecordData(BitMemoryWriter& metaWriter_,
										 BitMemoryWriter& dnaWriter_,
										 BitMemoryWriter& quaWriter_,
										 BitMemoryWriter& headWriter_,
										 const BinPackSettings& settings_,
										 BinaryBinDescriptor& binDesc_,
										 const FastqRecord& rec_)
{
	if (!settings_.hasConstLen)
		metaWriter_.PutBits(rec_.seqLen - settings_.minLen, settings_.bitsPerLen);

	StoreNextRecord(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, rec_);

	//ASSERT((uint64)binDesc_.rawDnaSize + (uint64)rec_.seqLen < (1ULL << 32));
	//ASSERT((uint64)binDesc_.rawHeadSize + (uint64)rec_.headLen < (1ULL << 32));

	binDesc_.rawDnaSize += rec_.seqLen;
	binDesc_.rawHeadSize += rec_.headLen;
}


void FastqNodesPackerSE::ReadRecordData(BitMemoryReader& metaReader_,
										BitMemoryReader& dnaReader_,
										BitMemoryReader& quaReader_,
										BitMemoryReader& headReader_,
										const BinPackSettings& settings_,
										FastqChunk& fqChunk_,
										FastqRecord& rec_)
{
	// read the sequence length
	//
	if (settings_.hasConstLen)
		rec_.seqLen = settings_.minLen;
	else
		rec_.seqLen = metaReader_.GetBits(settings_.bitsPerLen) + settings_.minLen;
	ASSERT(rec_.seqLen > 0 && rec_.seqLen < FastqRecord::MaxSeqLen);;


	// bind the sequence to data chunk
	//
	rec_.seq = (char*)(fqChunk_.data.Pointer() + fqChunk_.size);
	rec_.qua = rec_.seq + rec_.seqLen;

	if (binConfig.archiveType.readsHaveHeaders)
	{
		rec_.head = rec_.qua + rec_.seqLen;
	}


	// read the sequence
	//
	bool b = ReadNextRecord(metaReader_, dnaReader_, quaReader_, headReader_, settings_, rec_);
	ASSERT(b);
	ASSERT(fqChunk_.size + (uint64)(rec_.seqLen * 2 + rec_.headLen) <= fqChunk_.data.Size());

	if (settings_.suffixLen > 0)
	{
		std::copy(settings_.signatureString.data(),
				  settings_.signatureString.data() + binConfig.minimizer.signatureLen,
				  rec_.seq + rec_.minimPos);
	}

	// update the chunk size
	//
	fqChunk_.size += (uint64)(rec_.seqLen * 2 + rec_.headLen);			// * 2 <-- for quality
}


void FastqNodesPackerSE::StoreExactMatch(BitMemoryWriter& metaWriter_,
										 BitMemoryWriter& /*dnaWriter_*/,
										 BitMemoryWriter& quaWriter_,
										 BitMemoryWriter& headWriter_,
										 const BinPackSettings& settings_,
										 BinaryBinDescriptor& binDesc_,
										 const FastqRecord& rec_)
{
	metaWriter_.PutBit(rec_.IsReadReverse());

	StoreQuality(metaWriter_, quaWriter_, settings_, rec_);

	if (binConfig.archiveType.readsHaveHeaders)
	{
		StoreHeader(metaWriter_, headWriter_, settings_, rec_);
	}

	//ASSERT((uint64)binDesc_.rawDnaSize + (uint64)rec_.seqLen < (1ULL << 32));
	//ASSERT((uint64)binDesc_.rawHeadSize + (uint64)rec_.headLen < (1ULL << 32));

	binDesc_.rawDnaSize += rec_.seqLen;
	binDesc_.rawHeadSize += rec_.headLen;
}


void FastqNodesPackerSE::ReadExactMatch(BitMemoryReader& metaReader_,
										BitMemoryReader& /*dnaReader_*/,
										BitMemoryReader& quaReader_,
										BitMemoryReader& headReader_,
										const BinPackSettings& settings_,
										FastqChunk& fqChunk_,
										const FastqRecord& mainRec_,
										FastqRecord& rec_)
{
	bool isRev = metaReader_.GetBit() != 0;
	rec_.SetReadReverse(isRev);

	rec_.seq = (char*)(fqChunk_.data.Pointer() + fqChunk_.size);
	rec_.seqLen = mainRec_.seqLen;
	rec_.qua = rec_.seq + rec_.seqLen;


	// copy the record contents from the paren
	//
	std::copy(mainRec_.seq, mainRec_.seq + mainRec_.seqLen, rec_.seq);
	rec_.minimPos = mainRec_.minimPos;

	ReadQuality(metaReader_, quaReader_, settings_, rec_);

	if (binConfig.archiveType.readsHaveHeaders)
	{
		rec_.head = rec_.qua + rec_.seqLen;
		ReadHeader(metaReader_, headReader_, settings_, rec_);
	}


	fqChunk_.size += (uint64)(rec_.seqLen * 2 + rec_.headLen);
}


// PE packer
//
void FastqNodesPackerPE::StoreRecordData(BitMemoryWriter& metaWriter_,
										 BitMemoryWriter& dnaWriter_,
										 BitMemoryWriter& quaWriter_,
										 BitMemoryWriter& headWriter_,
										 const BinPackSettings& settings_,
										 BinaryBinDescriptor& binDesc_,
										 const FastqRecord& rec_)
{
	if (!settings_.hasConstLen)
	{
		metaWriter_.PutBits(rec_.seqLen - settings_.minLen, settings_.bitsPerLen);
		metaWriter_.PutBits(rec_.auxLen - settings_.minLen, settings_.bitsPerLen);
	}

	if (settings_.suffixLen != 0)
		metaWriter_.PutBit(rec_.IsPairSwapped());

	FastqRecord pair = rec_.GetPair();
	StoreNextRecord(metaWriter_, dnaWriter_, quaWriter_, headWriter_, settings_, rec_);
	StoreNextRecord(metaWriter_, dnaWriter_, quaWriter_, headWriter_, defaultPairSettings, pair);

	//ASSERT((uint64)binDesc_.rawDnaSize + (uint64)rec_.seqLen + pair.seqLen < (1ULL << 32));
	//ASSERT((uint64)binDesc_.rawHeadSize + (uint64)rec_.headLen < (1ULL << 32));

	binDesc_.rawDnaSize += rec_.seqLen + pair.seqLen;
	binDesc_.rawHeadSize += rec_.headLen;
}


void FastqNodesPackerPE::ReadRecordData(BitMemoryReader& metaReader_,
										BitMemoryReader& dnaReader_,
										BitMemoryReader& quaReader_,
										BitMemoryReader& headReader_,
										const BinPackSettings& settings_,
										FastqChunk& fqChunk_,
										FastqRecord& rec_)
{
	// read the sequence length
	//
	if (settings_.hasConstLen)
	{
		rec_.seqLen = settings_.minLen;
		rec_.auxLen = rec_.seqLen;
	}
	else
	{
		rec_.seqLen = metaReader_.GetBits(settings_.bitsPerLen) + settings_.minLen;
		rec_.auxLen = metaReader_.GetBits(settings_.bitsPerLen) + settings_.minLen;
	}
	ASSERT(rec_.seqLen > 0 && rec_.seqLen < FastqRecord::MaxSeqLen);

	if (settings_.suffixLen)
		rec_.SetPairSwapped(metaReader_.GetBit() != 0);


	// bind the sequence to data chunk
	//
	rec_.seq = (char*)(fqChunk_.data.Pointer() + fqChunk_.size);
	rec_.qua = rec_.seq + rec_.seqLen + rec_.auxLen;

	if (binConfig.archiveType.readsHaveHeaders)
	{
		rec_.head = rec_.qua + rec_.seqLen + rec_.auxLen;
	}


	// read the _1
	//
	bool b = ReadNextRecord(metaReader_, dnaReader_, quaReader_, headReader_, settings_, rec_);
	ASSERT(b);
	ASSERT(fqChunk_.size + (uint64)((rec_.seqLen + rec_.auxLen) * 2 + rec_.headLen) <= fqChunk_.data.Size());

	if (settings_.suffixLen > 0)
	{
		std::copy(settings_.signatureString.data(),
				  settings_.signatureString.data() + binConfig.minimizer.signatureLen,
				  rec_.seq + rec_.minimPos);
	}

	// read the _2
	//
	FastqRecord pair = rec_.GetPair();
	b = ReadNextRecord(metaReader_, dnaReader_, quaReader_, headReader_, defaultPairSettings, pair);
	ASSERT(b);

	// update the chunk size
	//
	fqChunk_.size += (uint64)((rec_.seqLen + pair.seqLen) * 2  + rec_.headLen);		// * 2  <-- for quality
}


void FastqNodesPackerPE::StoreExactMatch(BitMemoryWriter& metaWriter_,
										 BitMemoryWriter& dnaWriter_,
										 BitMemoryWriter& quaWriter_,
										 BitMemoryWriter& headWriter_,
										 const BinPackSettings& settings_,
										 BinaryBinDescriptor& binDesc_,
										 const FastqRecord& rec_)
{
	metaWriter_.PutBit(rec_.IsReadReverse());
	metaWriter_.PutBit(rec_.IsPairSwapped());

	StoreQuality(metaWriter_, quaWriter_, settings_, rec_);

	if (binConfig.archiveType.readsHaveHeaders)
	{
		StoreHeader(metaWriter_, headWriter_, settings_, rec_);
	}

	FastqRecord pair = rec_.GetPair();
	StoreNextRecord(metaWriter_, dnaWriter_, quaWriter_, headWriter_, defaultPairSettings, pair);

	//ASSERT((uint64)binDesc_.rawDnaSize + (uint64)(rec_.seqLen + pair.seqLen) < (1ULL << 32));
	//ASSERT((uint64)binDesc_.rawHeadSize + (uint64)rec_.headLen < (1ULL << 32));

	binDesc_.rawDnaSize += rec_.seqLen + pair.seqLen;
	binDesc_.rawHeadSize += rec_.headLen;
}


void FastqNodesPackerPE::ReadExactMatch(BitMemoryReader& metaReader_,
										BitMemoryReader& dnaReader_,
										BitMemoryReader& quaReader_,
										BitMemoryReader& headReader_,
										const BinPackSettings& settings_,
										FastqChunk& fqChunk_,
										const FastqRecord& mainRec_,
										FastqRecord& rec_)
{
	// read SE
	//
	bool isRev = metaReader_.GetBit() != 0;
	bool isSwap = metaReader_.GetBit() != 0;
	rec_.SetReadReverse(isRev);
	rec_.SetPairSwapped(isSwap);

	rec_.seq = (char*)(fqChunk_.data.Pointer() + fqChunk_.size);
	rec_.seqLen = mainRec_.seqLen;
	rec_.auxLen = mainRec_.auxLen;
	rec_.qua = rec_.seq + rec_.seqLen + rec_.auxLen;

	// copy the parent sequence and read own quality
	//
	std::copy(mainRec_.seq, mainRec_.seq + mainRec_.seqLen, rec_.seq);
	rec_.minimPos = mainRec_.minimPos;

	ReadQuality(metaReader_, quaReader_, settings_, rec_);

	if (binConfig.archiveType.readsHaveHeaders)
	{
		rec_.head = rec_.qua + rec_.seqLen + rec_.auxLen;
		ReadHeader(metaReader_, headReader_, settings_, rec_);
	}


	// read PE
	//
	FastqRecord pair = rec_.GetPair();
	bool b = ReadNextRecord(metaReader_, dnaReader_, quaReader_, headReader_, defaultPairSettings, pair);
	ASSERT(b);


	// update the chunk size
	//
	fqChunk_.size += (uint64)((rec_.seqLen + pair.seqLen) * 2 + rec_.headLen);
}




// dynamic memory
//
void IFastqNodesPackerDyn::UnpackFromBin(const BinaryBinBlock& binBin_,
									 std::vector<FastqRecord>& reads_,
									 GraphEncodingContext& graph_,
									 FastqRecordBinStats& stats_,
									 IFastqChunkCollection& fqChunk_,
									 bool append_)
{
	ASSERT(binBin_.blockType == BinaryBinBlock::SingleSignatureType);
	ASSERT(binBin_.auxDescriptors.size() > 0);
	ASSERT(binBin_.signature != 0);

	ASSERT(!append_);

	// calculate global bin properties
	//
	uint64 recordsCount = 0;

	for (const auto& desc : binBin_.auxDescriptors)
	{
		ASSERT(desc.recordsCount > 0);
		recordsCount += (uint64)desc.recordsCount;
	}

	ASSERT(recordsCount < (1 << 28));
	ASSERT(recordsCount != 0);

	// initialize buffers
	//
	uint64 recId = 0;

	// reset stats for initialization
	//
	ASSERT(reads_.size() == 0);

	stats_.Clear();

	std::vector<uint64> chunkSizes;
	uint64 curSize = 0;
	for (const BinaryBinDescriptor& desc : binBin_.auxDescriptors)
	{
		uint64 sliceSize = desc.rawHeadSize + desc.rawDnaSize * 2;

		if (curSize + sliceSize > DefaultMaxChunkSize)
		{
			chunkSizes.push_back(curSize + sliceSize);
			curSize = 0;
		}
		else
		{
			curSize += sliceSize;
		}
	}
	if (curSize > 0)
	{
		chunkSizes.push_back(curSize);
		curSize = 0;
	}

	auto chunkSizesIter = chunkSizes.cbegin();


	/*
	for (DataChunk* dc : fqChunk_.chunks)
	{
		if (dc->data.Size() < DefaultMaxChunkSize)
			dc->data.Extend(DefaultMaxChunkSize);
		dc->size = 0;
	}
	*/

	if (fqChunk_.chunks.size() != 0)
		fqChunk_.Clear();

	//ASSERT(fqChunk_.chunks.size() == 0);

	fqChunk_.chunks.push_back(new DataChunk(*chunkSizesIter++));
	DataChunk* currentChunk = fqChunk_.chunks.back();

	reads_.resize(recordsCount);


#if EXTRA_MEM_OPT
	reads_.shrink_to_fit();
#endif


	// initialize IO
	//
	BitMemoryReader metaReader(binBin_.metaData, binBin_.metaSize);
	BitMemoryReader dnaReader(binBin_.dnaData, binBin_.dnaSize);
	BitMemoryReader quaReader(binBin_.quaData, binBin_.quaSize);
	BitMemoryReader headReader(binBin_.headData, binBin_.headSize);

	// initialize settings
	//
	BinPackSettings settings;
	settings.signatureId = binBin_.signature;

	if (binBin_.signature != binConfig.minimizer.TotalMinimizersCount())	// N bin?
	{
		settings.suffixLen = binConfig.minimizer.signatureLen;
		binConfig.minimizer.GenerateMinimizer(binBin_.signature, settings.signatureString.data());
	}
	else
	{
		settings.suffixLen = 0;
	}
	settings.usesHeaders = binConfig.archiveType.readsHaveHeaders;


	// unpack records
	//
	for (const BinaryBinDescriptor& desc : binBin_.auxDescriptors)
	{
		uint64 sliceSize = desc.rawHeadSize + desc.rawDnaSize * 2;

		// todo: optimize chunk selection
		if (currentChunk->size + sliceSize > currentChunk->data.Size())
		{
			ASSERT(chunkSizesIter != chunkSizes.end());
			fqChunk_.chunks.push_back(new DataChunk(*chunkSizesIter++));
			currentChunk = fqChunk_.chunks.back();
		}


		ASSERT(desc.recordsCount > 0);

		const uint64 initialMetaPos = metaReader.Position();
		const uint64 initialDnaPos = dnaReader.Position();
		const uint64 initialQuaPos = quaReader.Position();
		const uint64 initialHeadPos = headReader.Position();

		// read bin header
		//
		settings.minLen = metaReader.GetBits(LenBits);
		settings.maxLen = metaReader.GetBits(LenBits);
		ASSERT(settings.minLen > 0);
		ASSERT(settings.maxLen >= settings.minLen);

		settings.hasReadGroups = metaReader.GetBit() != 0;
		settings.hasConstLen = (settings.minLen == settings.maxLen);

		if (!settings.hasConstLen)
		{
			settings.bitsPerLen = bit_length(settings.maxLen - settings.minLen);
			ASSERT(settings.bitsPerLen > 0);
		}

		// read records from bin
		//
		uint64 totalReads = recId + (uint64)desc.recordsCount;

		const uint64 initialChunkSize = currentChunk->size;
		while (recId < totalReads)
		{
			// todo: store the number of total nodes to avoid copying with resizing
			graph_.nodes.push_back(MatchNode());
			MatchNode& node = graph_.nodes.back();

			ReadNextNode(graph_, node,  reads_,  recId,
						 metaReader, dnaReader, quaReader, headReader,
						 settings, *currentChunk);
		}

		metaReader.FlushInputWordBuffer();
		dnaReader.FlushInputWordBuffer();
		quaReader.FlushInputWordBuffer();
		headReader.FlushInputWordBuffer();

		ASSERT(metaReader.Position() - initialMetaPos == (uint64)desc.metaSize);
		ASSERT(dnaReader.Position() - initialDnaPos == (uint64)desc.dnaSize);
		ASSERT(quaReader.Position() - initialQuaPos == (uint64)desc.quaSize);
		ASSERT(headReader.Position() - initialHeadPos == (uint64)desc.headSize);

		ASSERT(currentChunk->size == initialChunkSize + (uint64)desc.rawDnaSize * 2 + (uint64)desc.rawHeadSize);
	}


	// update the stats
	//
	stats_.minSeqLen = settings.minLen;
	stats_.maxSeqLen = settings.maxLen;

	/*
#if EXTRA_MEM_OPT

	std::vector<DataChunk*> validChunks;
	for (DataChunk* dc : fqChunk_.chunks)
	{
		if (dc->size > 0)
			validChunks.push_back(dc);
		else
			delete dc;
	}
	std::swap(fqChunk_.chunks, validChunks);

#endif
	*/

}
