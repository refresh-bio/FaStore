/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "../core/Globals.h"

#include <iostream>
#include <string.h>

#include "main.h"
#include "CompressorModule.h"
#include "Params.h"

#include "../core/Utils.h"
#include "../core/Thread.h"
#include "../core/version.h"


uint32 InputArguments::AvailableCoresNumber = mt::thread::hardware_concurrency();
uint32 InputArguments::DefaultThreadNumber = MIN(8, InputArguments::AvailableCoresNumber);


int main(int argc_, const char* argv_[])
{
	if (argc_ < 1 + 3 || (argv_[1][0] != 'e' && argv_[1][0] != 'd'))
	{
		usage();
		return -1;
	}

	InputArguments args;
	if (!parse_arguments(argc_, argv_, args))
		return -1;

	if (args.mode == InputArguments::EncodeMode)
		return bin2dnarch(args);
	return dnarch2dna(args);
}


void usage()
{
	std::cerr << "\n\n\t\t--- FaStore ---\n\n\n";
	std::cerr << "fastore_pack -- FASTQ reads compression tool\n";
	std::cerr << "Version: " << GetAppVersion() << " @ (" << GetCompilationTime() << ")\n";
	std::cerr << "Authors:  Lukasz Roguski\n          Idoia Ochoa\n          Mikel Hernaez\n          Sebastian Deorowicz\n\n\n";

	std::cerr << "usage:\tfastore_pack <e|d> [options] -i<input_file> -o<output_file>\n";

	std::cerr << "\nI/O options:\n";
	std::cerr << "\t-i<file>\t: input file(s) prefix";
	std::cerr << "\t-o<file>\t: output files prefix\n";
	std::cerr << "\t-o\"<f1> <f2> ... <fn>\": output FASTQ files list (PE mode)" << '\n';
	std::cerr << "\t-z\t\t: use paired-end mode, default: false\n";

	std::cerr << "\nrecords LZ-matching options:\n";
	std::cerr << "\t-f<n>\t\t: minimum bin size to filter, default: " << BinExtractorParams::Default::MinBinSize << '\n';
	std::cerr << "\t-e<n>\t\t: encode threshold value, default: 0 (auto)\n";
	std::cerr << "\t-m<n>\t\t: mismatch cost, default: " << ReadsClassifierParams::Default::MismatchCost << '\n';
	std::cerr << "\t-s<n>\t\t: shift cost, default: " << ReadsClassifierParams::Default::ShiftCost << '\n';
	std::cerr << "\t-w<n>\t\t: max LZ match window, default: " << ReadsClassifierParams::Default::MaxLzWindowSize << '\n';

	std::cerr << "\t-r\t\t: reduce Hard Reads by extra search in prefix buffer, default: false \n";
	std::cerr << "\t-l\t\t: reduce Expensive LZ-matches by extra search in prefix buffer, default: false \n";

	std::cerr << "\nrecords LZ-matching options in paired-end mode:\n";
	std::cerr << "\t-E<n>\t\t: pair encode threshold value, default: 0 (auto)\n";
	std::cerr << "\t-W<n>\t\t: max LZ match window, default: " << ReadsClassifierParams::Default::MaxPairLzWindowSize << '\n';

	std::cerr << "\nmatch tree and consensus building options:\n";
	std::cerr << "\t-c<n>\t\t: min consensus size, default: " << ReadsContigBuilderParams::Default::MinConsensusSize << '\n';
	std::cerr << "\t-q<n>\t\t: max record shift, default: 0 (auto)\n";
	std::cerr << "\t-n<n>\t\t: max new variants per read, default: " << ReadsContigBuilderParams::Default::MaxNewVariantsPerRead << '\n';
	std::cerr << "\t-d<n>\t\t: max Hamming distance, default: " << ReadsContigBuilderParams::Default::MaxHammingDistance << '\n';

	//std::cerr << "\t-x\t\t: use stored sub-tree topology while compressing, default: false\n";
    
    std::cerr << "QVZ Options are:\n\n";
    std::cerr << "\t-T\t\t: Target average distortion, measured as specified by -d or -D (default 1)\n";
	std::cerr << "\t-D <M|L|A>\t: Optimize for MSE, Log(1+L1), L1 distortions, respectively (default: MSE)\n";
	std::cerr << "\t-M<FILE>\t: Optimize using the custom distortion matrix specified in FILE\n";
	std::cerr << "\t-U<FILE>\t: Write the uncompressed lossy values to FILE (default: off)\n";
	std::cerr << "\t-F\t\t: output full FASTQ reads instead only q-scores (used with -U), default: false\n";

    std::cerr << "\nFor custom distortion matrices, a 72x72 matrix of values must be provided as the cost of reconstructing\n";
    std::cerr << "the x-th row as the y-th column, where x and y range from 0 to 71 (inclusive) corresponding to the possible Phred scores.\n";


	std::cerr << "\ngeneral options:\n";
	std::cerr << "\t-t<n>\t\t: threads count, default: " << InputArguments::DefaultThreadNumber << '\n';
	std::cerr << "\t-v\t\t: verbose mode, default: false\n";
}


int bin2dnarch(const InputArguments& args_)
{
	try
	{
		if (args_.pairedEndMode)
		{
			CompressorModulePE module;
			module.Bin2Dnarch(args_.inputFile, args_.outputFiles[0], args_.params, args_.auxParams,
							  args_.threadsNum, args_.verboseMode);
		}
		else
		{
			CompressorModuleSE module;
			module.Bin2Dnarch(args_.inputFile, args_.outputFiles[0], args_.params, args_.auxParams,
					args_.threadsNum, args_.verboseMode);
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "Error: " << e.what() << std::endl;
		return -1;
	}

	return 0;
}


int dnarch2dna(const InputArguments& args_)
{
	try
	{
		if (args_.pairedEndMode)
		{
			CompressorModulePE module;
			module.Dnarch2Dna(args_.inputFile, args_.outputFiles[0], args_.outputFiles[1], args_.threadsNum);
		}
		else
		{
			CompressorModuleSE module;
			module.Dnarch2Dna(args_.inputFile, args_.outputFiles[0], args_.threadsNum);
		}

		// very robust way to oputput re-shuffled reads to tmp output file
		//
		if (args_.params.quality.method == QualityCompressionParams::MET_QVZ && args_.auxParams.dry_run)
		{
			ASSERT(args_.auxParams.f_uncompressed != NULL);
			fclose(args_.auxParams.f_uncompressed);

			if (args_.auxParams.pe_mutex != NULL)
				delete args_.auxParams.pe_mutex;

			if (args_.auxParams.f_uncompressed_2 != NULL)
				fclose(args_.auxParams.f_uncompressed_2);
		}
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
	outArgs_.mode = (argv_[1][0] == 'e') ? InputArguments::EncodeMode : InputArguments::DecodeMode;

    
    // DEFAULT QVZ OPTIONS
    outArgs_.qvzOpts.verbose = 0;
    outArgs_.qvzOpts.stats = 0;
    outArgs_.qvzOpts.distortion = DISTORTION_MSE;
    outArgs_.qvzOpts.D = 1;
    outArgs_.qvzOpts.uncompressed_name = NULL;
    
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
		{
			pval = to_num((const uchar*)param + 2, len - 2);
			ASSERT(pval >= 0);
		}
		const char* str = param + 2;
		const uint32 slen = len - 2;

		switch (param[1])
		{
			case 'i':	outArgs_.inputFile.assign(str, str + slen);		break;
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

			case 't':	outArgs_.threadsNum = pval;									break;
			case 'v':
            {
                outArgs_.verboseMode = true;
                outArgs_.qvzOpts.stats = 1;
                outArgs_.qvzOpts.verbose = 1;
                break;
                
            }
			case 'z':	outArgs_.pairedEndMode = true;								break;

			case 'f':	outArgs_.params.extractor.minBinSize = pval;				break;

			case 'w':	outArgs_.params.classifier.maxLzWindowSize = pval;			break;
			case 'W':	outArgs_.params.classifier.maxPairLzWindowSize = pval;		break;
			case 'e':	outArgs_.params.classifier.encodeThresholdValue = pval;		break;
			case 'E':	outArgs_.params.classifier.pairEncodeThresholdValue = pval;	break;
			case 's':	outArgs_.params.classifier.shiftCost = pval;				break;
			case 'm':	outArgs_.params.classifier.mismatchCost = pval;				break;

			case 'r':	outArgs_.params.classifier.extraReduceHardReads = true;		break;
			case 'l':	outArgs_.params.classifier.extraReduceExpensiveLzMatches = true; break;

			//case 'x':	outArgs_.params.useStoredTopology = true;					break;

			case 'q':	outArgs_.params.consensus.maxRecordShiftDifference = pval;	break;
			case 'n':	outArgs_.params.consensus.maxNewVariantsPerRead = pval;		break;
			case 'd':	outArgs_.params.consensus.maxHammingDistance = pval;		break;
			case 'c':	outArgs_.params.consensus.minConsensusSize = pval;			break;

			// QVZ
            case 'U':
            {
				outArgs_.auxParams.dry_run = true;
                outArgs_.qvzOpts.uncompressed = 1;
                //outArgs_.params.quality.uncompressed_name = &(param[2]);

				int beg = 2;
				for (int i = 2; i < len-1; ++i)
				{
					if (param[i] == ' ' || param[i] == '\n')
					{
						outArgs_.auxParams.uncompressed_filename = std::string(param + beg, i - beg);
						beg = i + 1;
					}
				}
				if (len - beg > 0)
				{
					outArgs_.auxParams.uncompressed_filename_2 = std::string(param + beg, len - beg);

					if (outArgs_.auxParams.uncompressed_filename.length() == 0)
						std::swap(outArgs_.auxParams.uncompressed_filename, outArgs_.auxParams.uncompressed_filename_2);
				}

                break;
            }
			case 'F':	outArgs_.auxParams.output_fastq = true;						break;

            case 'M':
            {
                outArgs_.qvzOpts.distortion = DISTORTION_CUSTOM;
                outArgs_.qvzOpts.dist_file = &(param[2]);
                break;
            }
            case 'T':   outArgs_.qvzOpts.D = atof(&(param[2]));                      break;
            case 'D':
                switch (param[2]) {
                    case 'M':
                        outArgs_.qvzOpts.distortion = DISTORTION_MSE;
                        break;
                    case 'L':
                        outArgs_.qvzOpts.distortion = DISTORTION_LORENTZ;
                        break;
                    case 'A':
                        outArgs_.qvzOpts.distortion = DISTORTION_MANHATTAN;
                        break;
                    default:
                        printf("Prebuilt distortion measure not supported, using MSE.\n");
                        break;
                }
                break;
                
                
                

		}
	}

	// check params
	//
	if (outArgs_.inputFile.size() == 0)
	{
		std::cerr << "Error: no input file specified\n";
		return false;
	}

	if (outArgs_.outputFiles.size() == 0)
	{
		std::cerr << "Error: no output file(s) specified\n";
		return false;
	}

	if (outArgs_.mode == InputArguments::DecodeMode)
	{
		if (outArgs_.pairedEndMode && outArgs_.outputFiles.size() != 2)
		{
			std::cerr << "Error: no output file(s) specified for PE mode\n";
			return false;
		}
	}

	if (outArgs_.threadsNum == 0 || outArgs_.threadsNum > 64)
	{
		std::cerr << "Error: invalid number of threads specified\n";
		return false;
	}

	if (outArgs_.auxParams.dry_run)
	{
		if (outArgs_.auxParams.uncompressed_filename.length() == 0)
		{
			std::cerr << "Error: no output file specified for dry run mode\n";
			return false;
		}

		if (outArgs_.pairedEndMode)
		{
			if (outArgs_.auxParams.uncompressed_filename_2.length() == 0)
			{
				std::cerr << "Error: no output file #2 specified for dry run mode\n";
				return false;
			}

			if((outArgs_.auxParams.f_uncompressed_2 = fopen(outArgs_.auxParams.uncompressed_filename_2.c_str(), "w")) == NULL)
			{
				std::cerr << "ERROR: output file #2 for dry run couldn't be opened"<< std::endl;
				return false;
			}

			outArgs_.auxParams.pe_mutex = new std::mutex();
		}


		if((outArgs_.auxParams.f_uncompressed = fopen(outArgs_.auxParams.uncompressed_filename.c_str(), "w")) == NULL)
		{
			std::cerr << "ERROR: output file #1 for dry run couldn't be opened"<< std::endl;
			return false;
		}
	}

	return true;
}
