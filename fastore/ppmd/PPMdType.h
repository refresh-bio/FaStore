/****************************************************************************
 *  This file is part of PPMd library                                       *
 *  Original code designed and distributed to public domain                 *
 *  by Dmitry Shkarin 1997, 1999-2001, 2006, 2013                           *
 *  Contents: compilation parameters and miscelaneous definitions           *
 *  Comments: system & compiler dependent file                              *
 ****************************************************************************/

#ifndef _PPMDTYPE_H_
#define _PPMDTYPE_H_
//#define NDEBUG

#include <assert.h>

#define PPMD_ERR	-1
#define PPMD_OK		0

//#define _WIN32_ENVIRONMENT_
//#define _DOS32_ENVIRONMENT_
#define _POSIX_ENVIRONMENT_
//#define _UNKNOWN_ENVIRONMENT_
#if defined(_WIN32_ENVIRONMENT_)+defined(_DOS32_ENVIRONMENT_)+defined(_POSIX_ENVIRONMENT_)+defined(_UNKNOWN_ENVIRONMENT_) != 1
#error Only one environment must be defined
#endif /* defined(_WIN32_ENVIRONMENT_)+defined(_DOS32_ENVIRONMENT_)+defined(_POSIX_ENVIRONMENT_)+defined(_UNKNOWN_ENVIRONMENT_) != 1 */

#if defined(_WIN32_ENVIRONMENT_)
#include <windows.h>
#else /* _DOS32_ENVIRONMENT_ || _POSIX_ENVIRONMENT_ || _UNKNOWN_ENVIRONMENT_ */
typedef int   BOOL;
#define FALSE 0
#define TRUE  1
typedef unsigned char  BYTE;                // it must be equal to uint8_t
typedef unsigned short WORD;                // it must be equal to uint16_t
typedef unsigned int  DWORD;                // it must be equal to uint32_t
typedef unsigned long QWORD;                // it must be equal to uint64_t
typedef unsigned int   UINT;
#endif /* defined(_WIN32_ENVIRONMENT_)  */
                        /* Optimal definitions for processors:              */
//#define _32_NORMAL    /* IA-32                                            */
#define _64_NORMAL    /* AMD64/EM64T                                      */
//#define _32_EXOTIC    /* with request for 32bit alignment for uint32_t    */
//#define _64_EXOTIC    /* some unknown to me processors                    */
#if defined(_32_NORMAL)+defined(_64_NORMAL)+defined(_32_EXOTIC)+defined(_64_EXOTIC) != 1
#error Only one processor type must be defined
#endif /* defined(_32_NORMAL)+defined(_64_NORMAL)+defined(_32_EXOTIC)+defined(_64_EXOTIC) != 1 */

#define _PAD_TO_64(Dummy)
#if defined(_32_NORMAL) || defined(_64_NORMAL)
typedef BYTE  _BYTE;
typedef WORD  _WORD;
typedef DWORD _DWORD;
#else
#pragma message ("Warning: real memory usage will be twice larger")
typedef WORD  _BYTE;
typedef DWORD _WORD;
#if defined(_32_EXOTIC)
typedef DWORD _DWORD;
#undef _PAD_TO_64
#define _PAD_TO_64(Dummy) DWORD Dummy;
#else
typedef QWORD _DWORD;
#endif
#endif /* defined(_32_NORMAL) || defined(_64_NORMAL) */

#if !defined(NDEBUG)
BOOL TestCompilation();                     /* for testing our data types   */
#endif /* !defined(NDEBUG) */

#if !defined(_UNKNOWN_ENVIRONMENT_) && !defined(__GNUC__)
#define _FASTCALL __fastcall
#define _STDCALL  __stdcall
#else
#define _FASTCALL
#define _STDCALL
#endif /* !defined(_UNKNOWN_ENVIRONMENT_) && !defined(__GNUC__) */

/* _USE_THREAD_KEYWORD macro must be defined at compilation for creation    *
 * of multithreading applications. Some compilers generate correct code     *
 * with the use of standard '__thread' keyword (GNU C), some others use it  *
 * in non-standard way (BorlandC) and some use __declspec(thread) keyword   *
 * (IntelC, VisualC).                                                       */
#define _USE_THREAD_KEYWORD

#if defined(_USE_THREAD_KEYWORD)
#if defined(_MSC_VER)
#define _THREAD
#define _THREAD1 __declspec(thread)
#elif defined(__GNUC__)
#define _THREAD
#define _THREAD1 __thread
#else /* __BORLANDC__ */
#define _THREAD __thread
#define _THREAD1
#endif /* defined(_MSC_VER) */
#else
#define _THREAD
#define _THREAD1
#endif /* defined(_USE_THREAD_KEYWORD) */

const DWORD PPMdSignature=0x84ACAF8F;
enum { PROG_VAR='J', MAX_O=16 };            /* maximum allowed model order  */
#define _USE_PREFETCHING                    /* it gives 2-6% speed gain     */

template <class T>
inline T CLAMP(const T& X,const T& LoX,const T& HiX) { return (X >= LoX)?((X <= HiX)?(X):(HiX)):(LoX); }
template <class T>
inline void SWAP(T& t1,T& t2) { T tmp=t1; t1=t2; t2=tmp; }

#if 0
/* PPMd module works with file streams via ...GETC/...PUTC macros only      */
typedef FILE _PPMD_FILE;
#define _PPMD_E_GETC(fp)   getc(fp)
#define _PPMD_E_PUTC(c,fp) putc((c),fp)
#define _PPMD_D_GETC(fp)   getc(fp)
#define _PPMD_D_PUTC(c,fp) putc((c),fp)
#else
class ByteStream;
typedef ByteStream _PPMD_FILE;
#define _PPMD_E_GETC(bs)   bs->Get()
#define _PPMD_E_PUTC(c,bs) bs->Put(c)
#define _PPMD_D_GETC(bs)   bs->Get()
#define _PPMD_D_PUTC(c,bs) bs->Put(c)
#endif

#endif /* !defined(_PPMDTYPE_H_) */
