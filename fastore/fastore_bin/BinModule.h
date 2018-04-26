/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_BINMODULE
#define H_BINMODULE

#include "../core/Globals.h"

#include <string>
#include <vector>

#include "Params.h"


/**
 * A standalone modules for binning/un-binning single/paired-end FASTQ data
 *
 */
class BinModuleSE
{
public:
	void Fastq2Bin(const std::vector<std::string>& inFastqFiles_,
				   const std::string& outBinFile_,
				   const BinModuleConfig& config_,
				   uint32 threadNum_ = 1,
				   bool compressedInput_ = false, bool verboseMode_ = false);

	void Bin2Dna(const std::string& inBinFile_,
				 const std::string& outFile_);
};


class BinModulePE
{
public:
	void Fastq2Bin(const std::vector<std::string>& inFastqFiles_1_,
				   const std::vector<std::string>& inFastqFiles_2_,
				   const std::string& outBinFile_,
				   const BinModuleConfig& config_,
				   uint32 threadNum_ = 1, bool compressedInput_ = false, bool verboseMode_ = false);

	void Bin2Dna(const std::string& inBinFile_,
				 const std::string& outFile_1_,
				 const std::string& outFile_2_);
};


#endif // H_BINMODULE
