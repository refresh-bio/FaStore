/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_MAIN
#define H_MAIN

#include "../core/Globals.h"

#include <string>
#include <vector>

#include "Params.h"


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

	BinBalanceParameters params;

	bool useMatePairs;
	uint32 threadsNum;
	bool verboseMode;

	std::vector<std::string> inputFiles;
	std::vector<std::string> outputFiles;

	InputArguments()
		:	mode(EncodeMode)
		,	useMatePairs(false)
		,	threadsNum(DefaultThreadNumber)
		,	verboseMode(DefaultVerboseMode)
	{}
};


void usage();
int bin2bin(const InputArguments& args_);
int bin2dna(const InputArguments& args_);
bool parse_arguments(int argc_, const char* argv_[], InputArguments& outArgs_);

int main(int argc_, const char* argv_[]);


#endif //H_ MAIN
