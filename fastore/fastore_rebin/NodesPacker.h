/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef NODESPACKER_H
#define NODESPACKER_H


#include "../fastore_bin/Globals.h"
#include "../fastore_bin/FastqPacker.h"
#include "../fastore_bin/Node.h"

#include <deque>


/**
 * Binary-encodes the connection graph between matched reads
 * and the reads themselves
 */
class IFastqNodesPacker : public IFastqPacker
{
public:
	using IFastqPacker::IFastqPacker;

	// these methods operate on DnaRecord* data pointers, where the actual DnaRecord content is stored
	// elsewhere in the memory (used after DNA records distribution into bins)
	void PackToBins(const std::map<uint32, MatchNodesPtrBin>& dnaBins_,
					BinaryBinBlock& binBins_);


	// TODO: we need to re-check use-cases of this one
	void PackToBin(const GraphEncodingContext& graph_,
				   BinaryBinBlock& binBin_,
				   uint32 signatureId_);


	void UnpackFromBin(const BinaryBinBlock& binBin_,
					   std::vector<FastqRecord>& reads_,
					   GraphEncodingContext& graph_,
					   FastqRecordBinStats& stats_,
					   FastqChunk& fqChunk_,
						   bool append_ = true);

protected:
		virtual void StoreRecordData(BitMemoryWriter& metaWriter_,
									 BitMemoryWriter& dnaWriter_,
									 BitMemoryWriter& quaWriter_,
									 BitMemoryWriter& headWriter_,
									 const BinPackSettings& settings_,
									 BinaryBinDescriptor& binDesc_,
									 const FastqRecord& rec_) = 0;

		virtual void ReadRecordData(BitMemoryReader& metaReader_,
									BitMemoryReader& dnaReader_,
									BitMemoryReader& quaReader_,
									BitMemoryReader& headReader_,
									const BinPackSettings& settings_,
									FastqChunk& fqChunk_,
									FastqRecord& rec_) = 0;

		virtual void StoreExactMatch(BitMemoryWriter& metaWriter_,
									 BitMemoryWriter& dnaWriter_,
									 BitMemoryWriter& quaWriter_,
									 BitMemoryWriter& headWriter_,
									 const BinPackSettings& settings_,
									 BinaryBinDescriptor& binDesc_,
									 const FastqRecord& rec_) = 0;

		virtual void ReadExactMatch(BitMemoryReader& metaReader_,
									BitMemoryReader& dnaReader_,
									BitMemoryReader& quaReader_,
									BitMemoryReader& headReader_,
									const BinPackSettings& settings_,
									FastqChunk& fqChunk_,
									const FastqRecord& mainRec_,
									FastqRecord& rec_) = 0;

protected:
		void PackToBin(const MatchNodesPtrBin& dnaBin_,
					   BitMemoryWriter& metaWriter_,
					   BitMemoryWriter& dnaWriter_,
					   BitMemoryWriter& quaWriter_,
					   BitMemoryWriter& headWriter_,
					   BinaryBinDescriptor& binDesc_,
					   bool nBin_ = false);

		void StoreNextNode(BitMemoryWriter& metaWriter_,
						   BitMemoryWriter& dnaWriter_,
						   BitMemoryWriter& quaWriter_,
						   BitMemoryWriter& headWriter_,
						   const BinPackSettings& settings_,
						   BinaryBinDescriptor& binDesc_,
						   const MatchNode* node_,
						   std::deque<const MatchNode*>& packContext_);

		void StoreNextGroup(BitMemoryWriter& metaWriter_,
							BitMemoryWriter& dnaWriter_,
							BitMemoryWriter& quaWriter_,
							BitMemoryWriter& headWriter_,
							const BinPackSettings& settings_,
							BinaryBinDescriptor& binDesc_,
							uint32 signatureId_,
							uint32 signaturePos_,
							uint64 groupSize_,
							std::deque<const MatchNode*>& group_);

		void ReadNextNode(GraphEncodingContext& graph_,
						  MatchNode& curNode_,
						  std::vector<FastqRecord>& reads_,
						  uint64& recIdx_,
						  BitMemoryReader& metaReader_,
						  BitMemoryReader& dnaReader_,
						  BitMemoryReader& quaReader_,
						  BitMemoryReader& headReader_,
						  const BinPackSettings& settings_,
						  FastqChunk& fqChunk_);
};


class FastqNodesPackerSE : public virtual IFastqNodesPacker
{
public:
	//using IFastqNodesPacker::IFastqNodesPacker;
	FastqNodesPackerSE(const BinModuleConfig& binConfig_)
		:	IFastqNodesPacker(binConfig_)
	{}

protected:
	void StoreRecordData(BitMemoryWriter& metaWriter_,
								 BitMemoryWriter& dnaWriter_,
								 BitMemoryWriter& quaWriter_,
								 BitMemoryWriter& headWriter_,
								 const BinPackSettings& settings_,
								 BinaryBinDescriptor& binDesc_,
								 const FastqRecord& rec_);

	void ReadRecordData(BitMemoryReader& metaReader_,
								BitMemoryReader& dnaReader_,
								BitMemoryReader& quaReader_,
								BitMemoryReader& headReader_,
								const BinPackSettings& settings_,
								FastqChunk& fqChunk_,
								FastqRecord& rec_);

	void StoreExactMatch(BitMemoryWriter& metaWriter_,
								 BitMemoryWriter& dnaWriter_,
								 BitMemoryWriter& quaWriter_,
								 BitMemoryWriter& headWriter_,
								 const BinPackSettings& settings_,
								 BinaryBinDescriptor& binDesc_,
								 const FastqRecord& rec_);

	void ReadExactMatch(BitMemoryReader& metaReader_,
								BitMemoryReader& dnaReader_,
								BitMemoryReader& quaReader_,
								BitMemoryReader& headReader_,
								const BinPackSettings& settings_,
								FastqChunk& fqChunk_,
								const FastqRecord& mainRec_,
								FastqRecord& rec_);
};



class FastqNodesPackerPE : public virtual IFastqNodesPacker
{
public:
		FastqNodesPackerPE(const BinModuleConfig& binConfig_)
				:	IFastqNodesPacker(binConfig_)
		{
			defaultPairSettings.minLen = defaultPairSettings.maxLen = 1;		// some default dummy value
			defaultPairSettings.hasConstLen = true;
			defaultPairSettings.usesHeaders = false;
		}

protected:
		void StoreRecordData(BitMemoryWriter& metaWriter_,
									 BitMemoryWriter& dnaWriter_,
									 BitMemoryWriter& quaWriter_,
									 BitMemoryWriter& headWriter_,
									 const BinPackSettings& settings_,
									 BinaryBinDescriptor& binDesc_,
									 const FastqRecord& rec_);

		void ReadRecordData(BitMemoryReader& metaReader_,
									BitMemoryReader& dnaReader_,
									BitMemoryReader& quaReader_,
									BitMemoryReader& headReader_,
									const BinPackSettings& settings_,
									FastqChunk& fqChunk_,
									FastqRecord& rec_);

		void StoreExactMatch(BitMemoryWriter& metaWriter_,
									 BitMemoryWriter& dnaWriter_,
									 BitMemoryWriter& quaWriter_,
									 BitMemoryWriter& headWriter_,
									 const BinPackSettings& settings_,
									 BinaryBinDescriptor& binDesc_,
									 const FastqRecord& rec_);

		void ReadExactMatch(BitMemoryReader& metaReader_,
									BitMemoryReader& dnaReader_,
									BitMemoryReader& quaReader_,
									BitMemoryReader& headReader_,
									const BinPackSettings& settings_,
									FastqChunk& fqChunk_,
									const FastqRecord& mainRec_,
									FastqRecord& rec_);
private:
		BinPackSettings defaultPairSettings;
};



// dynamic memory
//
class IFastqNodesPackerDyn : public virtual IFastqNodesPacker
{
public:
	//using IFastqNodesPacker::IFastqNodesPacker;
	IFastqNodesPackerDyn(const BinModuleConfig& binConfig_)
		:	IFastqNodesPacker(binConfig_)
	{}

	void UnpackFromBin(const BinaryBinBlock& binBin_,
					   std::vector<FastqRecord>& reads_,
					   GraphEncodingContext& graph_,
					   FastqRecordBinStats& stats_,
					   IFastqChunkCollection& fqChunk_,
					   bool append_ = true);

protected:
	static const uint64 DefaultMaxChunkSize = 32 << 20;
};


class FastqNodesPackerDynSE : public FastqNodesPackerSE
							, public virtual IFastqNodesPackerDyn
{
public:
	//using FastqNodesPackerSE::FastqNodesPackerSE;
	FastqNodesPackerDynSE(const BinModuleConfig& binConfig_)
		:	IFastqNodesPacker(binConfig_)
		,	IFastqNodesPackerDyn(binConfig_)
		,	FastqNodesPackerSE(binConfig_)
	{}
};


class FastqNodesPackerDynPE : public FastqNodesPackerPE
							, public virtual IFastqNodesPackerDyn
{
public:
	//using FastqNodesPackerPE::FastqNodesPackerPE;
	FastqNodesPackerDynPE(const BinModuleConfig& binConfig_)
		:	IFastqNodesPacker(binConfig_)
		,	IFastqNodesPackerDyn(binConfig_)
		,	FastqNodesPackerPE(binConfig_)
	{}
};


#endif // FASTQNODESPACKER_H
