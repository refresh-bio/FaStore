/****************************************************************************
 *  This file is part of PPMd library                                       *
 *  Contents: PPMd library wrapper                                          *
 *  Author: Lucas Roguski, 2014                                             *
 ****************************************************************************/

#ifndef _PPMD_H_
#define _PPMD_H_

#include <stdint.h>

#ifdef  __cplusplus
extern "C" {
#endif

class PpmdEncoder
{
public:
	static const unsigned int DefaultMemorySizeMb = 32;
	static const unsigned int DefaultOrder = 4;

	PpmdEncoder();
	~PpmdEncoder();

	static bool Encode(unsigned char *inBuffer_, uint64_t inBufferSize_,
					   unsigned char* outBuffer_, uint64_t &outBufferSize_,
					   unsigned int order_ = DefaultOrder, unsigned int memorySizeMb_ = DefaultMemorySizeMb);

	bool StartCompress(unsigned int order_ = DefaultOrder, unsigned int memorySize_ = DefaultMemorySizeMb);
	bool EncodeNextMember(unsigned char* inBuffer_, uint64_t inBufferSize_,
						  unsigned char* outBuffer_, uint64_t &outBufferSize_);
	bool FinishCompress();

private:
	unsigned int order;
};

class PpmdDecoder
{
public:
	static const unsigned int DefaultMemorySizeMb = PpmdEncoder::DefaultMemorySizeMb;

	PpmdDecoder();
	~PpmdDecoder();

	static bool Decode(unsigned char* inBuffer_, uint64_t inBufferSize_,
					   unsigned char* outBuffer_, uint64_t &outBufferSize_,
					   unsigned int memorySizeMb_ = DefaultMemorySizeMb);

	bool StartDecompress(unsigned int memorySizeMb_ = DefaultMemorySizeMb);
	bool DecodeNextMember(unsigned char* inBuffer_, uint64_t inBufferSize_,
						  unsigned char* outBuffer_, uint64_t &outBufferSize_);
	bool FinishDecompress();
};

#ifdef  __cplusplus
}
#endif

#endif // _PPMD_H_
