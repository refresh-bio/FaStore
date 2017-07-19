/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_MAIN
#define H_MAIN

#include "../fastore_bin/Globals.h"

#include <string>
#include <vector>

#include "Params.h"
#include "codebook.h"



struct InputArguments
{
	enum ModeEnum
	{
		EncodeMode,
		DecodeMode
	};

	static const bool DefaultVerboseMode = false;

	static uint32 AvailableCoresNumber;
	static uint32 DefaultThreadNumber;

	ModeEnum mode;

	std::string inputFile;
	std::vector<std::string> outputFiles;

	CompressorParams params;
	CompressorAuxParams auxParams;

	uint32 threadsNum;
	bool verboseMode;
	bool pairedEndMode;
    
    // QVZ options
    struct qv_options_t qvzOpts;
    
	InputArguments()
		:	threadsNum(DefaultThreadNumber)
		,	verboseMode(DefaultVerboseMode)
		,	pairedEndMode(false)
	{}
};


void usage();
int bin2dnarch(const InputArguments& args_);
int dnarch2dna(const InputArguments& args_);
bool parse_arguments(int argc_, const char* argv_[], InputArguments& outArgs_);

int main(int argc_, const char* argv_[]);


#endif // H_MAIN
