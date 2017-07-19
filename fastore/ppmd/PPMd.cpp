/****************************************************************************
 *  This file is part of PPMd library                                       *
 *  Contents: PPMd library wrapper                                          *
 *  Author: Lucas Roguski, 2014                                             *
 ****************************************************************************/

#include "PPMd.h"

#include <stdlib.h>
#include <cstring>

#include "PPMdType.h"
#include "Stream.hpp"

#define PPMD_CONTROL_BYTE			0xCA
#define PPMD_DEFAULT_ORDER			PpmdEncoder::DefaultOrder
#define PPMD_DEFAULT_ALLOC_SIZE_MB	PpmdEncoder::DefaultMemorySizeMb

// Declarations for C-wrappers:
//
int ppmd_start_suballocator(unsigned int subAllocatorSize_);
int ppmd_stop_suballocator();
unsigned int ppmd_used_memory();

int ppmd_compress(unsigned char* pInMemory_, uint64_t inSize_,
				  unsigned char* pOutMemory_, uint64_t* outSize_,
				  unsigned int maxOrder_ = PPMD_DEFAULT_ORDER, bool useAutoAllocation_ = true,
				  unsigned int allocatorSizeMb_ = PPMD_DEFAULT_ALLOC_SIZE_MB, bool doOrderCutOff_ = false);

int ppmd_decompress(unsigned char* pInMemory_, uint64_t inSize_,
					unsigned char* pOutMemory_, uint64_t* outSize_,
					bool useAutoAllocation_ = true, unsigned int allocatorSizeMb_ = PPMD_DEFAULT_ALLOC_SIZE_MB);


// C++ class wrapper:
//
PpmdEncoder::PpmdEncoder()
	:	order(0)
{}

PpmdEncoder::~PpmdEncoder()
{}

bool PpmdEncoder::Encode(unsigned char *inBuffer_, uint64_t inBufferSize_,
						 unsigned char *outBuffer_, uint64_t &outBufferSize_,
						 unsigned int order_, unsigned int memorySizeMb_)
{
	return ppmd_compress(inBuffer_, inBufferSize_,
						 outBuffer_, &outBufferSize_,
						 order_, true, memorySizeMb_) == PPMD_OK;
}

bool PpmdEncoder::StartCompress(unsigned  order_, unsigned int memorySize_)
{
	order = order_;
	return ppmd_start_suballocator(memorySize_) == PPMD_OK;
}

bool PpmdEncoder::EncodeNextMember(unsigned char *inBuffer_, uint64_t inBufferSize_,
								   unsigned char *outBuffer_, uint64_t &outBufferSize_)
{
	return ppmd_compress(inBuffer_, inBufferSize_,
						 outBuffer_, &outBufferSize_,
						 order, false) == PPMD_OK;
}

bool PpmdEncoder::FinishCompress()
{
	return ppmd_stop_suballocator() == PPMD_OK;
}


PpmdDecoder::PpmdDecoder()
{}

PpmdDecoder::~PpmdDecoder()
{}

bool PpmdDecoder::Decode(unsigned char *inBuffer_, uint64_t inBufferSize_,
						 unsigned char *outBuffer_, uint64_t &outBufferSize_,
						 unsigned int memorySizeMb_)
{
	return ppmd_decompress(inBuffer_, inBufferSize_,
						   outBuffer_, &outBufferSize_,
						   true, memorySizeMb_) == PPMD_OK;
}

bool PpmdDecoder::StartDecompress(unsigned int memorySizeMb_)
{
	return ppmd_start_suballocator(memorySizeMb_) == PPMD_OK;
}

bool PpmdDecoder::DecodeNextMember(unsigned char *inBuffer_, uint64_t inBufferSize_,
								   unsigned char *outBuffer_, uint64_t &outBufferSize_)
{
	return ppmd_decompress(inBuffer_, inBufferSize_,
						   outBuffer_, &outBufferSize_,
						   false) == PPMD_OK;
}

bool PpmdDecoder::FinishDecompress()
{
	return ppmd_stop_suballocator() == PPMD_OK;
}


// C-style function wrappers:
//

// SubAlloc:
UINT _STDCALL GetUsedMemory();
int _STDCALL StopSubAllocator();
int _STDCALL StartSubAllocator(UINT SASize);

// Model:
void EncodeFile(_PPMD_FILE* EncodedFile,_PPMD_FILE* DecodedFile, int MaxOrder,BOOL CutOff);
void DecodeFile(_PPMD_FILE* DecodedFile,_PPMD_FILE* EncodedFile, int MaxOrder,BOOL CutOff);

int ppmd_compress(unsigned char *pInMemory_, uint64_t inSize_,
				  unsigned char *pOutMemory_, uint64_t* outSize_,
				  unsigned int maxOrder_, bool useAutoAllocation_,
				  unsigned int allocatorSizeMb_, bool doOrderCutOff_)
{
	if (maxOrder_ == 0 || maxOrder_ > 16
			|| outSize_ == NULL
			|| *outSize_ == 0 || *outSize_ < inSize_
			|| pInMemory_ == NULL || pOutMemory_ == NULL)
		return PPMD_ERR;

	ByteStream streamIn(pInMemory_, inSize_);
	ByteStream streamOut(pOutMemory_, *outSize_);

	// write only important header shit
	//
	unsigned int header = ((int)doOrderCutOff_ << 7) | maxOrder_;
	streamOut.Put(PPMD_CONTROL_BYTE);
	streamOut.Put(header);

	if (useAutoAllocation_)
	{
		ppmd_start_suballocator(allocatorSizeMb_);
	}

	EncodeFile(&streamOut, &streamIn, maxOrder_, doOrderCutOff_);

	*outSize_ = streamOut.Position();

	if (useAutoAllocation_)
	{
		ppmd_stop_suballocator();
	}

	return PPMD_OK;
}

int ppmd_decompress(unsigned char *pInMemory_, uint64_t inSize_,
					unsigned char *pOutMemory_, uint64_t* outSize_,
					bool useAutoAllocation_, unsigned int allocatorSizeMb_)
{
	if (outSize_ == 0 || *outSize_ < inSize_
			|| pInMemory_ == NULL || pOutMemory_ == NULL)
		return PPMD_ERR;

	ByteStream streamIn(pInMemory_, inSize_);
	ByteStream streamOut(pOutMemory_, *outSize_);

	// read the header shit
	//
	unsigned int header = streamIn.Get();
	if (header != PPMD_CONTROL_BYTE)
		return PPMD_ERR;

	header = streamIn.Get();
	bool doCoutoff = (header & 80) != 0;
	int maxOrder = header & 0x7F;

	if (maxOrder == 0 || maxOrder > 16)
		return PPMD_ERR;

	if (useAutoAllocation_)
	{
		ppmd_start_suballocator(allocatorSizeMb_);	// default: 10 MB
	}

	DecodeFile(&streamOut, &streamIn, maxOrder, doCoutoff);

	*outSize_ = streamOut.Position();

	if (useAutoAllocation_)
	{
		ppmd_stop_suballocator();
	}

	return PPMD_OK;
}

int ppmd_start_suballocator(unsigned int subAllocatorSize_)
{
	return StartSubAllocator(subAllocatorSize_);
}

int ppmd_stop_suballocator()
{
	return StopSubAllocator();
}

unsigned int ppmd_used_memory()
{
	return GetUsedMemory();
}
