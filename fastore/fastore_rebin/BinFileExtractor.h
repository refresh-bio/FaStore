/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINFILEEXTRACTOR
#define H_BINFILEEXTRACTOR

#include "../core/Globals.h"
#include "../core/BinBlockData.h"

#include "../fastore_bin/BinFile.h"


struct BinExtractorParams
{
	struct Default
	{
		static const uint32 MinBinSize = 256;				// 64 --> 512
	};

	uint32 minBinSize;

	BinExtractorParams()
		:	minBinSize(Default::MinBinSize)
	{}
};


/**
 * Extracts the bins by prrovided signature id
 *
 */
class BinFileExtractor : public BinFileReader
{
public:
	static const uint32 DefaultMinimumBinSize = 64;

	BinFileExtractor(uint32 minBinSize_ = DefaultMinimumBinSize);

	void StartDecompress(const std::string& fileName_, BinModuleConfig& params_);

	bool ExtractNextSmallBin(BinaryBinBlock& bin_);
	bool ExtractNextStdBin(BinaryBinBlock& bin_);
	bool ExtractNBin(BinaryBinBlock& bin_);

	using BinFileReader::ReadBlock;

	uint64 BlockCount() const
	{
		return fileHeader.blockCount;
	}

	std::map<uint32, const BinFileFooter::BinInfo*> GetBlockDescriptors(bool stdBlocks_ = true) const;
	std::pair<uint32, const BinFileFooter::BinInfo*> GetNBlockDescriptor() const;

private:
	using BinFileReader::ReadNextBlock;

	const uint32 minBinSize;

	std::vector<uint32> stdSignatures;
	std::vector<uint32> smallSignatures;
	std::vector<uint32>::const_iterator stdSignatureIterator;
	std::vector<uint32>::const_iterator smallSignatureIterator;
};


#endif // H_BINFILEEXTRACTOR
