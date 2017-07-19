/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_GLOBALS
#define H_GLOBALS


// debugging
//
#ifndef NDEBUG
#	define DEBUG 1
#endif

#define EXTRA_MEM_OPT 1

#if defined(DEBUG)
#	define DEV_DEBUG_MODE 0
#endif

// Visual Studio warning supression
//
#if defined (_WIN32)

#	define _CRT_SECURE_NO_WARNINGS
#	pragma warning(disable : 4996) // D_SCL_SECURE
#	pragma warning(disable : 4244) // conversion uint64 to uint32
#	pragma warning(disable : 4267)
#	pragma warning(disable : 4800) // conversion byte to bool
#endif


// assertions
//
#if defined(DEBUG) || defined(_DEBUG)
#   include "stdio.h"
#	include "assert.h"
#	define ASSERT(x) assert(x)
#   define SOFT_ASSERT(expr) \
	((expr)	\
	? __ASSERT_VOID_CAST (0) \
	: __ASSERT_VOID_CAST(fprintf(stderr, "Soft assertion failed: '%s' in '%s' @ %d : '%s'\n", __STRING(expr), __FILE__, __LINE__, __ASSERT_FUNCTION)))
#else
#	define ASSERT(x)
#endif

#define STATIC_ASSERT(X) static_assert((X), "static assertion fail")



// global constants
//
#define MAX_SIGNATURE_LEN 16


#define MAX_LZ_SE 255		// FIXME: a very bad workaround !!!
#define MAX_LZ_PE 4096


// global macros
//
#define BIT(x)							(1 << (x))
#define MIN(x,y)						((x) <= (y) ? (x) : (y))
#define MAX(x,y)						((x) >= (y) ? (x) : (y))
#define ABS(x)							((x) >=  0  ? (x) : (-x))
#define TFREE(ptr)						{ if (ptr != NULL) delete ptr; ptr = NULL; }


// basic types
//
#include <stdint.h>
typedef int8_t				int8;
typedef uint8_t				uchar, byte, uint8;
typedef int16_t				int16;
typedef uint16_t			uint16;
typedef int32_t				int32;
typedef uint32_t			uint32;
typedef int64_t				int64;
typedef uint64_t			uint64;




/**
 * A common operator interface used in multi-threaded mode
 *
 */
class IOperator
{
public:
	virtual ~IOperator() {}

	void operator() ()
	{
		Run();
	}

	virtual void Run() = 0;
};



/**
 * A common coder interface used for compression codecs
 *
 */
struct ICoder
{
	virtual ~ICoder() {}

	virtual void Start() = 0;
	virtual void End() = 0;
};



// TODO: split declarations depending on the project -- bin/pack
//
struct FastqRecord;
class BitMemoryReader;
class BitMemoryWriter;
class IFastqStreamReaderSE;
class IFastqStreamWriter;
class BinFileReader;
class BinFileWriter;
struct BinaryBinBlock;
struct DataChunk;
class Buffer;
struct MinimizerParameters;
struct CategorizerParameters;
struct BinModuleConfig;

struct GraphEncodingContext;
struct ExactMatchesGroup;
struct ContigDefinition;
struct TreeTransferDefinition;

struct FastqRecordBinStats;
struct FastqRawBlockStats;


#endif // H_GLOBALS
