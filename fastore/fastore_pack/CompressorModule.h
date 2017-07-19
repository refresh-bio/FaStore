/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_DNARCHMODULE
#define H_DNARCHMODULE

#include "../fastore_bin/Globals.h"


#include "../fastore_bin/FastqRecord.h"

#include <string>
#include <vector>

#include "Params.h"



/**
 * A standalone modules for compressing/de-compressing single/paired-end FASTQ data
 *
 */
class CompressorModuleSE
{
public:
	void Bin2Dnarch(const std::string& inBinFile_, const std::string& outArchiveFile_,
					const CompressorParams& params_, const CompressorAuxParams& auxParams_ = CompressorAuxParams(),
					uint32 threadsNum_ = 1, bool verboseMode_ = false);
	void Dnarch2Dna(const std::string& inArchiveFile_, const std::string& outDnaFile_,
					uint32 threadsNum_ = 1);
};


class CompressorModulePE
{
public:
	void Bin2Dnarch(const std::string& inBinFile_, const std::string& outArchiveFile_,
					const CompressorParams& params_, const CompressorAuxParams& auxParams_ = CompressorAuxParams(),
					uint32 threadsNum_ = 1, bool verboseMode_ = false);

	void Dnarch2Dna(const std::string& inArchiveFile_, const std::string& outDnaFile1_,
					const std::string& outDnaFile2_, uint32 threadsNum_ = 1);
};



#endif // H_DNARCHMODULE
