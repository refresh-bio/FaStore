/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DNAPACKER
#define H_DNAPACKER

#include "BinBlockData.h"
#include "FastqRecord.h"


/**
 * Binary-encodes binned FASTQ reads -- interface
 *
 */
class IFastqPacker
{
public:
	IFastqPacker(const BinModuleConfig& binConfig_);

	virtual ~IFastqPacker() {}

protected:
	struct BinPackSettings
	{
		static const uint32 MaxSignatureLength = 32;

		bool hasConstLen;
		bool hasReadGroups;			// TODO: use this for backward compatibility
		uint32 minLen;
		uint32 maxLen;
		uint32 suffixLen;
		uint32 bitsPerLen;			// TODO: variable length support

		// auxilary members
		uint32 signatureId;
		std::array<char, MaxSignatureLength> signatureString;
		bool usesHeaders;

		BinPackSettings()
			:	hasConstLen(false)
			,	hasReadGroups(false)
			,	minLen((uint32)-1)
			,	maxLen(0)
			,	suffixLen(0)
			,	bitsPerLen(0)
			,	signatureId(0)
			,	usesHeaders(false)
		{}
	};

	static const uint32 LenBits = 8;

	const BinModuleConfig& binConfig;

	std::array<char, 128> dnaToIdx;
	std::array<char, 8> idxToDna;

	std::array<char, 64> quaToIdx_8bin;
	std::array<char, 8> idxToQua_8bin;


	// TODO: initialize/finalize() to setup writers skipping passing by ref
	//


	void StoreNextRecord(BitMemoryWriter& metaWriter_,
						 BitMemoryWriter& dnaWriter_,
						 BitMemoryWriter& quaWriter_,
						 BitMemoryWriter& headWriter_,
						 const BinPackSettings& settings_,
						 const FastqRecord& rec_);

	// TODO: here we can also add computing of the FastqStats
	// to be able to obtain per-bin stats
	bool ReadNextRecord(BitMemoryReader& metaReader_,
						BitMemoryReader& dnaReader_,
						BitMemoryReader& quaReader_,
						BitMemoryReader& headReader_,
						const BinPackSettings& settings_,
						FastqRecord& rec_);

	void StoreDna(BitMemoryWriter& metaWriter_,
				  BitMemoryWriter& dnaWriter_,
				  const BinPackSettings& settings_,
				  const FastqRecord& rec_);
	void StoreQuality(BitMemoryWriter& metaWriter_,
					  BitMemoryWriter& quaWriter_,
					  const BinPackSettings& settings_,
					  const FastqRecord& rec_);
	void StoreHeader(BitMemoryWriter& metaWriter_,
					 BitMemoryWriter& headWriter_,
					 const BinPackSettings& settings_,
					 const FastqRecord& rec_);

	void ReadDna(BitMemoryReader& metaReader_,
				 BitMemoryReader& dnaReader_,
				 const BinPackSettings& settings_,
				 FastqRecord& rec_);
	void ReadQuality(BitMemoryReader& metaReader_,
					 BitMemoryReader& quaReader_,
					 const BinPackSettings& settings_,
					 FastqRecord& rec_);
	void ReadHeader(BitMemoryReader& metaReader_,
					 BitMemoryReader& headReader_,
					 const BinPackSettings& settings_,
					 FastqRecord& rec_);
};


// TODO: templatize
//
class FastqRecordsPackerSE: public IFastqPacker
{
public:
	using IFastqPacker::IFastqPacker;

	// these methods operate on DnaRecord* data pointers, where the actual DnaRecord content is stored
	// elsewhere in the memory (used after DNA records distribution into bins)
	void PackToBins(const std::map<uint32, FastqRecordsPtrBin>& dnaBins_,
					BinaryBinBlock& binBins_);

	void PackToBin(const std::vector<FastqRecord>& records_,
				   BinaryBinBlock& binBlock_,
				   uint32 binId_);

	void UnpackFromBin(const BinaryBinBlock& binBin_,
					   std::vector<FastqRecord>& reads_,
					   FastqChunk& fqChunk_,
					   bool append_ = false);

protected:
	void PackToBin(const FastqRecordsPtrBin& dnaBin_,
				   BitMemoryWriter& metaWriter_,
				   BitMemoryWriter& dnaWriter_,
				   BitMemoryWriter& quaWriter_,
				   BitMemoryWriter& headWriter_,
				   BinaryBinDescriptor& binDesc_,
				   bool nBin_ = false);

	virtual void StoreRecords(const std::vector<FastqRecord*>& records_,
							  const BinPackSettings& settings_,
							  BitMemoryWriter& metaWriter_,
							  BitMemoryWriter& dnaWriter_,
							  BitMemoryWriter& quaWriter_,
							  BitMemoryWriter& headWriter_,
							  BinaryBinDescriptor& binDesc_);

	virtual void ReadRecords(std::vector<FastqRecord>::iterator itBegin_,
							std::vector<FastqRecord>::iterator itEnd_,
							const BinPackSettings& settings_,
							BitMemoryReader& metaReader_,
							BitMemoryReader& dnaReader_,
							BitMemoryReader& quaReader_,
							BitMemoryReader& headReader_,
							FastqChunk& fqChunk_);
};


class FastqRecordsPackerPE: public FastqRecordsPackerSE
{
public:
	using FastqRecordsPackerSE::FastqRecordsPackerSE;

protected:
	void StoreRecords(const std::vector<FastqRecord*>& records_,
					  const BinPackSettings& settings_,
					  BitMemoryWriter& metaWriter_,
					  BitMemoryWriter& dnaWriter_,
					  BitMemoryWriter& quaWriter_,
					  BitMemoryWriter& headWriter_,
					  BinaryBinDescriptor& binDesc_);

	void ReadRecords(std::vector<FastqRecord>::iterator itBegin_,
					std::vector<FastqRecord>::iterator itEnd_,
					const BinPackSettings& settings_,
					BitMemoryReader& metaReader_,
					BitMemoryReader& dnaReader_,
					BitMemoryReader& quaReader_,
					BitMemoryReader& headReader_,
					FastqChunk& fqChunk_);
};


#endif // H_DNAPACKER
