/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../fastore_bin/Globals.h"

#include "BinFileExtractor.h"
#include "../fastore_bin/Exception.h"


BinFileExtractor::BinFileExtractor(uint32 minBinSize_)
	:	minBinSize(minBinSize_)
{}


void BinFileExtractor::StartDecompress(const std::string &fileName_, BinModuleConfig &params_)
{
	BinFileReader::StartDecompress(fileName_, params_);

	// skip the 'N' bin here
	//
	const auto nIter = fileFooter.binOffsets.find(fileFooter.params.minimizer.TotalMinimizersCount());

	for (auto iOff = fileFooter.binOffsets.begin(); iOff != nIter; iOff++)
	{
		if (iOff->second.totalRecordsCount >= minBinSize)
			stdSignatures.push_back(iOff->first);
		else
			smallSignatures.push_back(iOff->first);
	}

	// TODO: perform random shuffle of the descriptors
	//
	//std::random_shuffle(stdSignatures.begin(), stdSignatures.end());


	smallSignatureIterator = smallSignatures.begin();
	stdSignatureIterator = stdSignatures.begin();
}


bool BinFileExtractor::ExtractNextStdBin(BinaryBinBlock &bin_)
{
	bin_.Clear();
	if (stdSignatureIterator == stdSignatures.end())
		return false;

	ReadBlock(*stdSignatureIterator, &bin_);

	stdSignatureIterator++;
	return true;
}


bool BinFileExtractor::ExtractNextSmallBin(BinaryBinBlock &bin_)
{
	if (smallSignatureIterator == smallSignatures.end())
		return false;

	ReadBlock(*smallSignatureIterator, &bin_);

	smallSignatureIterator++;

	return true;
}


bool BinFileExtractor::ExtractNBin(BinaryBinBlock &bin_)
{
	const uint32 nBinId = fileFooter.params.minimizer.TotalMinimizersCount();
	if (fileFooter.binOffsets.count(nBinId) == 0)
		return false;

	bin_.Clear();
	ReadBlock(nBinId, &bin_);
	return true;
}


std::map<uint32, const IBinFile::BinFileFooter::BinInfo*> BinFileExtractor::GetBlockDescriptors(bool stdBlocks_) const
{
	std::map<uint32, const BinFileFooter::BinInfo*> descriptors;

	const std::vector<uint32>* sigs = (stdBlocks_) ? (&stdSignatures) : (&smallSignatures);
	for (uint32 sig : *sigs)
		descriptors.insert(std::make_pair(sig, &fileFooter.binOffsets.at(sig)));

	return descriptors;
}

std::pair<uint32, const IBinFile::BinFileFooter::BinInfo*> BinFileExtractor::GetNBlockDescriptor() const
{
	uint32 nBinId = fileFooter.params.minimizer.TotalMinimizersCount();
	if (fileFooter.binOffsets.count(nBinId) == 0)
		return std::make_pair(nBinId, (const BinFileFooter::BinInfo*)NULL);

	return std::make_pair(nBinId, &fileFooter.binOffsets.at(nBinId));
}
