/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_REBINMODULE
#define H_REBINMODULE

#include "../core/Globals.h"

#include <string>
#include <vector>

#include "../core/FastqRecord.h"

#include "Params.h"


/**
 * A standalone module for re-binning single/paired-end FASTQ data
 *
 */
class RebinModule
{
public:
	void Bin2Bin(const std::string& inBinFile_,
				 const std::string& outBinFile_,
				 const BinBalanceParameters& params_,
				 uint32 threadNum_ = 1,
				 bool verboseMode_ = false);

	void Bin2Dna(const std::string& inBinFile_,
				 const std::vector<std::string>& outFiles_);
};


#endif // H_BINMODULE
