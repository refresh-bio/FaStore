/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#include "Globals.h"

#include <iostream>
#include <string.h>

#include "main.h"
#include "BinModule.h"
#include "Utils.h"
#include "Thread.h"
#include "../version.h"
#include "QVZ.h"


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
		return fastq2bin(args);
	return bin2dna(args);
}


void usage()
{
	std::cerr << "\n\n\t\t--- FaStore ---\n\n\n";
	std::cerr << "fastore_bin -- FASTQ reads binning tool\n\n";
	std::cerr << "Version: " << GetAppVersion() << " @ (" << GetCompilationTime() << ")\n";
	std::cerr << "Authors:  Lukasz Roguski\n          Idoia Ochoa\n          Mikel Hernaez\n          Sebastian Deorowicz\n\n\n";

	std::cerr << "usage: \tfastore_bin <e|d> [options]\n";

	std::cerr << "single-end compression options:\n";
    std::cerr << "\t-i<f>\t: input file" << '\n';
	std::cerr << "\t-i\"<f1> [<f2> ...]\": input FASTQ files list" << '\n';
	std::cerr << "\t-o<f>\t\t: output file" << '\n';

	std::cerr << "paired-end compression options:\n";
	std::cerr << "\t-z\t\t: use paired-end mode, default: false\n";
	std::cerr << "\t-g\t\t: input compressed in .gz format\n";
	std::cerr << "\t-i\"<f1_1> [<f2_1> ...] <f1_2>] [<f2_2> ...]\": input FASTQ files list" << '\n';
	std::cerr << "\t-o\"<f1_1> <f2_1>\": output FASTQ files list (PE mode)" << '\n';

	std::cerr << "clustering options:\n";
	std::cerr << "\t-p<n>\t\t: signature length, default: " << MinimizerParameters::Default::SignatureLength << '\n';
	std::cerr << "\t-s<n>\t\t: skip-zone length, default: " << MinimizerParameters::Default::SkipZoneLength << '\n';
	//std::cerr << "\t-c<n>\t\t: signature cutoff mask bits, default: " << MinimizerParameters::Default::SignatureMaskCutoffBits << '\n';
	std::cerr << "\t-m<n>\t\t: mimimum block bin size, default: " << CategorizerParameters::DefaultMinimumPartialBinSize << '\n';

	std::cerr << "read identifiers compression options:\n";
	std::cerr << "\t-H\t\t: keep identifiers (see option: -C), default: false\n";
	std::cerr << "\t-C\t\t: skip comments (content after space), default: false\n";


	std::cerr << "quality processing options:\n";
	std::cerr << "\t-q<n>\t\t: quality compression method [0-3], default: " << QualityCompressionParams::Default::Method << '\n';
	std::cerr << "\t\t *0\t: lossless\n";
	std::cerr << "\t\t *1\t: binary thresholding (optional -w parameter)\n";
	std::cerr << "\t\t *2\t: Illumina 8 bins\n";
	std::cerr << "\t\t *3\t: QVZ\n";
	std::cerr << "\t-w<n>\t\t: quality compression threshold (see: -q1), default: " << (uint32)QualityCompressionParams::Default::MinBinaryFilterThreshold << '\n';
	std::cerr << "\t-I\t\t: use Phred+64 quality scale offset default: false (using Phred+33)\n";
    
    std::cerr << "QVZ Options are:\n\n";
    std::cerr << "\t-T\t\t: Target average distortion, measured as specified by -d or -D (default 1)\n";
	std::cerr << "\t-D <M|L|A>\t: Optimize for MSE, Log(1+L1), L1 distortions, respectively (default: MSE)\n";
	std::cerr << "\t-M<FILE>\t: Optimize using the custom distortion matrix specified in FILE\n";
	std::cerr << "\t-U<FILE>\t: Write the uncompressed lossy values to FILE (default: off)\n";
    
    std::cerr << "\nFor custom distortion matrices, a 72x72 matrix of values must be provided as the cost of reconstructing\n";
    std::cerr << "the x-th row as the y-th column, where x and y range from 0 to 71 (inclusive) corresponding to the possible Phred scores.\n";

	std::cerr << "performance options:\n";
	std::cerr << "\t-b<n>\t\t: FASTQ input buffer size (in MB), default: " << (BinModuleConfig::DefaultFastqBlockSize >> 20) << '\n';
	std::cerr << "\t-t<n>\t\t: worker threads number, default: " << InputArguments::DefaultThreadNumber << '\n';
	std::cerr << "\t-v\t\t: verbose mode, default: false\n";
}


int fastq2bin(const InputArguments& args_)
{
	try
	{
		if (args_.config.archiveType.readType == ArchiveType::READ_PE)
		{
			BinModulePE module;
			auto f1 = decltype(args_.inputFiles)(args_.inputFiles.begin(),
												 args_.inputFiles.begin() + args_.inputFiles.size()/2);

			auto f2 = decltype(args_.inputFiles)(args_.inputFiles.begin() + args_.inputFiles.size()/2,
												 args_.inputFiles.end());

			module.Fastq2Bin(f1, f2, args_.outputFiles[0], args_.config,
							 args_.threadsNum, args_.compressedInput, args_.verboseMode);
		}
		else
		{
			BinModuleSE module;
			module.Fastq2Bin(args_.inputFiles, args_.outputFiles[0],
							 args_.config, args_.threadsNum,
							 args_.compressedInput, args_.verboseMode);
		}
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
		// TODO: automatically deduce the archive type
		//

		if (args_.config.archiveType.readType == ArchiveType::READ_PE)
		{
			BinModulePE module;
			module.Bin2Dna(args_.inputFiles[0], args_.outputFiles[0], args_.outputFiles[1]);
		}
		else
		{
			BinModuleSE module;
			module.Bin2Dna(args_.inputFiles[0], args_.outputFiles[0]);
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
	if (argv_[1][0] != 'e' && argv_[1][0] != 'd')
	{
		std::cerr << "Error: invalid mode specified\n";
		return false;
	}
	outArgs_.mode = (argv_[1][0] == 'e') ? InputArguments::EncodeMode : InputArguments::DecodeMode;
	outArgs_.config.binningType = BinModuleConfig::BIN_RECORDS;
    
    // DEFAULT QVZ OPTIONS
    outArgs_.config.quaParams.qvzOpts.verbose = 0;
    outArgs_.config.quaParams.qvzOpts.stats = 0;
    outArgs_.config.quaParams.qvzOpts.distortion = DISTORTION_MSE;
    outArgs_.config.quaParams.qvzOpts.D = 1;


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
			// input
			//
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

			case 'g':	outArgs_.compressedInput = true;								break;
			case 'b':	outArgs_.config.fastqBlockSize = (uint64)pval << 20;			break;

			case 't':	outArgs_.threadsNum = pval;										break;
			case 'v':	outArgs_.verboseMode = true; outArgs_.config.quaParams.qvzOpts.stats = 1; outArgs_.config.quaParams.qvzOpts.verbose = 1;									break;

			case 'z':	outArgs_.config.archiveType.readType = ArchiveType::READ_PE;	break;

			// clustering
			//
			case 'p':	outArgs_.config.minimizer.signatureLen = pval;					break;
			case 's':	outArgs_.config.minimizer.skipZoneLen = pval;					break;
			//case 'c':	outArgs_.config.minimizer.signatureMaskCutoffBits = pval;		break;
			case 'm':	outArgs_.config.catParams.minBlockBinSize = pval;				break;


			// headers
			//
			case 'H':	outArgs_.config.archiveType.readsHaveHeaders = true;            break;
			case 'C':	outArgs_.config.headParams.preserveComments = false;			break;


			// quality
			//
			case 'q':	outArgs_.config.quaParams.method = pval;						break;
			case 'w':	outArgs_.config.quaParams.binaryThreshold = pval;				break;
			case 'I':	outArgs_.config.archiveType.qualityOffset = ArchiveType::Illumina64QualityOffset;				break;

            
			// QVZ
			//
            case 'M':
            {
                outArgs_.config.quaParams.qvzOpts.distortion = DISTORTION_CUSTOM;
                outArgs_.config.quaParams.qvzOpts.dist_file = &(param[2]);
                break;
            }
            case 'T':   outArgs_.config.quaParams.qvzOpts.D = atof(&(param[2]));                                      break;
            case 'D':
                switch (param[2]) {
                    case 'M':
                        outArgs_.config.quaParams.qvzOpts.distortion = DISTORTION_MSE;
                        break;
                    case 'L':
                        outArgs_.config.quaParams.qvzOpts.distortion = DISTORTION_LORENTZ;
                        break;
                    case 'A':
                        outArgs_.config.quaParams.qvzOpts.distortion = DISTORTION_MANHATTAN;
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
	if (outArgs_.inputFiles.size() == 0)
	{
		std::cerr << "Error: no input file(s) specified\n";
		return false;
	}

	if (outArgs_.outputFiles.size() == 0)
	{
		std::cerr << "Error: no output file specified\n";
		return false;
	}

	if (outArgs_.mode == InputArguments::EncodeMode)
	{
		if (outArgs_.config.archiveType.readType == ArchiveType::READ_PE)
		{
			if (outArgs_.inputFiles.size() % 2 != 0)
			{
				std::cerr << "Error: invalid number of input files specified in PE mode\n";
				return false;
			}
		}
	}
	else if (outArgs_.mode == InputArguments::DecodeMode)
	{
		if (outArgs_.config.archiveType.readType == ArchiveType::READ_PE)
		{
			if (outArgs_.outputFiles.size() % 2 != 0)
			{
				std::cerr << "Error: invalid number of output files specified in PE mode\n";
				return false;
			}
		}
	}

	if (outArgs_.threadsNum == 0 || outArgs_.threadsNum > 64)
	{
		std::cerr << "Error: invalid number of threads specified\n";
		return false;
	}

	return true;
}


