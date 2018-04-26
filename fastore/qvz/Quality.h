/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef QUALITY_H
#define QUALITY_H

#include "../core/Globals.h"
#include "../core/FastqParser.h"

/**
 * Options for the QVZ compression process
 */
struct qv_options_t {
    uint8_t verbose;
    uint8_t stats;
    uint8_t uncompressed;
    uint8_t distortion;
	const char *dist_file;				// shall we move those to std::string ???
	const char *uncompressed_name;		//
    double D;
};

struct QualityCompressionParams
{
	enum CompressionMethod
	{
		MET_NONE = 0,
		MET_BINARY,
		MET_8BIN,
		MET_QVZ
	};

	struct Default
	{
		static const CompressionMethod Method = MET_NONE;
		static const uint8 MinBinaryFilterThreshold = 20;

		static const uint8 MinThresholdValue = 6;
		static const uint8 MaxThresholdValue = 40;
	};

	byte method;
	uint8 binaryThreshold;		// optional

	// QVZ options -- TODO: elaborate on the subset of the data to be stored
    struct qv_options_t qvzOpts;

	QualityCompressionParams()
		:	method(Default::Method)
		,	binaryThreshold(Default::MinBinaryFilterThreshold)
	{}

	uint32 BitsPerBase() const
	{
		ASSERT(method < 4);

		uint32 bpb[] = {6, 1, 3, 6};	// 6 bits in RAW mode: max 64 (40+) values
		return bpb[method];
	}
};


#endif // QUALITY_H
