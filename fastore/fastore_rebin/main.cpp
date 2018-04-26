/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../core/Globals.h"

#include <iostream>
#include <string.h>

#include "main.h"
#include "RebinModule.h"

#include "../core/Utils.h"
#include "../core/Thread.h"
#include "../core/version.h"


uint32 InputArguments::AvailableCoresNumber = mt::thread::hardware_concurrency();
uint32 InputArguments::DefaultThreadNumber = MIN(8, InputArguments::AvailableCoresNumber);


int main(int argc_, const char* argv_[])
{
	if (argc_ < 1 + 2)
	{
		usage();
		return -1;
	}

	InputArguments args;
	if (!parse_arguments(argc_, argv_, args))
		return -1;

	if (args.mode == InputArguments::EncodeMode)
		return bin2bin(args);
	else
		return bin2dna(args);
}


void usage()
{
	std::cerr << "\n\n\t\t--- FaStore ---\n\n\n";
	std::cerr << "fastore_rebin -- FASTQ reads re-binning tool\n";
	std::cerr << "Version: " << GetAppVersion() << " @ (" << GetCompilationTime() << ")\n";
	std::cerr << "Authors:  Lukasz Roguski\n          Idoia Ochoa\n          Mikel Hernaez\n          Sebastian Deorowicz\n\n\n";

	std::cerr << "usage:\tfastore_rebin <e|d> [params]\n";
	std::cerr << "options:\n";

	std::cerr << "\t-i<file>\t: input bin" << '\n';
	std::cerr << "\t-o<f>\t\t: output bin" << '\n';

	std::cerr << "\t-z\t\t: use paired-end mode, default: false\n";

	std::cerr << "\nre-binning options:\n";
	std::cerr << "\t-p<n>\t\t: signature parity, default: " << BinBalanceParameters::Default::SignatureParity << '\n';
	std::cerr << "\t-x<n>\t\t: min bin size to extract, default: " << BinBalanceParameters::Default::MinBinSizeToExtract << '\n';
	std::cerr << "\t-y<n>\t\t: min bin size to categorize, default: " << BinBalanceParameters::Default::MinBinSizeToCategorize << '\n';
	std::cerr << "\t-q<n>\t\t: min tree size to store, default: " << BinBalanceParameters::Default::MinTreeSize << '\n';

	std::cerr << "\nrecords LZ-matching options:\n";
	std::cerr << "\t-e<n>\t\t: encode threshold value, default: 0 (auto)\n";
	std::cerr << "\t-m<n>\t\t: mismatch cost, default: " << BinBalanceParameters::Default::MismatchCost << '\n';
	std::cerr << "\t-s<n>\t\t: shift cost, default: " << BinBalanceParameters::Default::ShiftCost << '\n';
	std::cerr << "\t-w<n>\t\t: max LZ match window, default: " << BinBalanceParameters::Default::MaxLzWindowSize << '\n';

	std::cerr << "\t-r\t\t: reduce Hard Reads by extra search in prefix buffer, default: false \n";
	std::cerr << "\t-l\t\t: reduce Expensive LZ-matches by extra search in prefix buffer, default: false \n";

	std::cerr << "\ngeneral options:\n";
	std::cerr << "\t-t<n>\t\t: worker threads number, default: " << InputArguments::DefaultThreadNumber << '\n';
	std::cerr << "\t-v\t\t: verbose mode, default: false\n";
}


int bin2bin(const InputArguments& args_)
{
	try
	{
		RebinModule module;
		module.Bin2Bin(args_.inputFiles[0], args_.outputFiles[0],
			args_.params, args_.threadsNum, args_.verboseMode);
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}

	return 0;
}


int bin2dna(const InputArguments& args_)
{
	try
	{
		RebinModule module;
		module.Bin2Dna(args_.inputFiles[0], args_.outputFiles);

	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}

	return 0;
}


bool parse_arguments(int argc_, const char* argv_[], InputArguments& outArgs_)
{
	if (argv_[1][0] != 'e' && argv_[1][0] != 'd')
	{
		std::cerr << "Error: invalid mode specified\n";
		return false;
	}

	outArgs_.mode = (argv_[1][0] == 'e') ? InputArguments::EncodeMode : InputArguments::DecodeMode;

	// parse params
	//
	for (int i = 2; i < argc_; ++i)
	{
		const char* param = argv_[i];
		if (param[0] != '-')
			continue;

		int pval = -1;
		int len = strlen(param);
		if (len > 2 && len < 10)
			pval = to_num((const uchar*)param + 2, len - 2);

		switch (param[1])
		{
			case 'i':
			{
				int beg = 2;
				for (int i = 2; i < len-1; ++i)
				{
					if (param[i] == ' ' || param[i] == '\n')
					{
						outArgs_.inputFiles.push_back(std::string(param + beg, i - beg));
						beg = i + 1;
					}
				}
				if (len - beg > 0)
					outArgs_.inputFiles.push_back(std::string(param + beg, len - beg));

				break;
			}

			case 'o':
			{
				int beg = 2;
				for (int i = 2; i < len-1; ++i)
				{
					if (param[i] == ' ' || param[i] == '\n')
					{
						outArgs_.outputFiles.push_back(std::string(param + beg, i - beg));
						beg = i + 1;
					}
				}
				if (len - beg > 0)
					outArgs_.outputFiles.push_back(std::string(param + beg, len - beg));

				break;
			}

			case 'p':	outArgs_.params.signatureParity = (uint32)pval;				break;
			case 'x':	outArgs_.params.minBinSizeToExtract = (uint32)pval;			break;
			case 'y':	outArgs_.params.minBinSizeToCategorize = (uint32)pval;		break;
			case 'q':	outArgs_.params.minTreeSize = pval;							break;

			case 'e':	outArgs_.params.classifier.encodeThresholdValue = pval;		break;
			case 's':	outArgs_.params.classifier.shiftCost = pval;				break;
			case 'm':	outArgs_.params.classifier.mismatchCost = pval;				break;
			case 'w':	outArgs_.params.classifier.maxLzWindowSize = pval;			break;

			case 'r':	outArgs_.params.classifier.extraReduceHardReads = true;		break;
			case 'l':	outArgs_.params.classifier.extraReduceExpensiveLzMatches = true; break;

			case 't':	outArgs_.threadsNum = pval;									break;
			case 'v':	outArgs_.verboseMode = true;								break;
			case 'z':	outArgs_.useMatePairs = true;								break;
		}
	}

	// check params
	//
	if (outArgs_.inputFiles.size() == 0)
	{
		std::cerr << "Error: no input file(s) specified\n";
		return false;
	}

	if (outArgs_.outputFiles.size() == 0)
	{
		std::cerr << "Error: no output file(s) specified\n";
		return false;
	}

	if (outArgs_.threadsNum == 0 || outArgs_.threadsNum > 64)
	{
		std::cerr << "Error: invalid number of threads specified\n";
		return false;
	}

	return true;
}
