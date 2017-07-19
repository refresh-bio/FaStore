/****************************************************************************
 *  This file is part of PPMd library                                       *
 *  Contents: 'Carryless rangecoder' by Dmitry Subbotin                     *
 *  Comments: this implementation is claimed to be a public domain          *
 ****************************************************************************/

enum { TOP=1 << 24, BOT=1 << 15 };
static _THREAD1 struct SUBRANGE { DWORD low, high, scale; } _THREAD Range;
static _THREAD1 DWORD _THREAD low, _THREAD code, _THREAD range;

inline void rcInitEncoder() { low=0; range=DWORD(-1); }
#define RC_ENC_NORMALIZE(stream) {                                          \
    while ((low ^ (low+range)) < TOP || range < BOT &&                      \
            ((range= -low & (BOT-1)),1)) {                                  \
        _PPMD_E_PUTC(low >> 24,stream);                                     \
        range <<= 8;                        low <<= 8;                      \
    }                                                                       \
}
inline void rcEncodeSymbol()
{
    low += Range.low*(range/=Range.scale);  range *= Range.high-Range.low;
}
inline void rcFlushEncoder(_PPMD_FILE* stream)
{
    for (UINT i=0;i < 4;i++) {
        _PPMD_E_PUTC(low >> 24,stream);     low <<= 8;
    }
}
inline void rcInitDecoder(_PPMD_FILE* stream)
{
    low=code=0;                             range=DWORD(-1);
    for (UINT i=0;i < 4;i++)
            code=(code << 8) | _PPMD_D_GETC(stream);
}
#define RC_DEC_NORMALIZE(stream) {                                          \
    while ((low ^ (low+range)) < TOP || range < BOT &&                      \
            ((range= -low & (BOT-1)),1)) {                                  \
        code=(code << 8) | _PPMD_D_GETC(stream);                            \
        range <<= 8;                        low <<= 8;                      \
    }                                                                       \
}
inline UINT rcGetCurrentCount() { return (code-low)/(range /= Range.scale); }
inline void rcRemoveSubrange()
{
    low += range*Range.low;                 range *= Range.high-Range.low;
}
inline UINT rcBinStart(UINT f0,UINT Shift)  { return f0*(range >>= Shift); }
inline UINT rcBinDecode  (UINT tmp)         { return (code-low >= tmp); }
inline void rcBinCorrect0(UINT tmp)         { range=tmp; }
inline void rcBinCorrect1(UINT tmp,UINT f1) { low += tmp;   range *= f1; }
