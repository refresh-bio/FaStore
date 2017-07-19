/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_UTILS
#define H_UTILS

#include "Globals.h"

#include <vector>
#include <map>
#include <string>


template <uint32 _TBitNum>
struct TBitMask
{
	static const uint64 Value = ((uint64)1 << (_TBitNum - 1)) | TBitMask<_TBitNum-1>::Value;
};


template <>
struct TBitMask<0>
{
	static const uint64 Value = 0;
};


template <uint32 _TNum>
struct TIsPow2
{
	static const bool Value = (_TNum & 1) == 0 && TIsPow2< (_TNum >> 1) >::Value;
};


template <>
struct TIsPow2<1>
{
	static const bool Value = true;
};


template <uint64 _TNum, uint32 _TPow>
struct TPow
{
	static const uint64 Value = _TNum * TPow<_TNum, _TPow - 1>::Value;
};


template <uint64 _TNum>
struct TPow<_TNum, 0>
{
	static const uint64 Value = 1;
};


template <uint32 _TNum>
struct TLog2
{
	static const uint32 Value = TLog2< (_TNum >> 1) >::Value + 1;
};


template <>
struct TLog2<1>
{
	static const uint32 Value = 0;
};
template <> struct TLog2<0> {};


template <uint32 _TMaxLen>
uint32 TStrLen(const char* str_)
{
	for (uint32 i = 0; i < _TMaxLen; ++i)
	{
		if (str_[i] == '\0')
			return i;
	}
	return _TMaxLen;
}


template <typename _T>
inline void TFree(_T*& p_)
{
	if (p_ != (_T*)0)
	{
		delete p_;
		p_ = (_T*)0;
	}
}


inline uint32 to_string(char* str, uint64 value)
{
	uint32 digits;
	uint64 power = 1;

	if (value == 0)
	{
		str[0] = '0';
		return 1;
	}

	for (digits = 0; digits < 16; ++digits)
	{
		if (value < power)
			break;
		power *= 10;
	}

	power /= 10;
	for (uint32 i = 0; power; ++i, power /= 10)
	{
		int32 d = value / power;
		str[i] = (uchar)('0' + d);
		value -= d * power;
	}

	return digits;
}

inline std::string to_string(uint32 value)
{
	if (value == 0)
		return "0";

	std::string str;
	uint32 digits;
	uint32 power = 1;
	for (digits = 0; digits < 10; ++digits)
	{
		if (value < power)
			break;
		power *= 10;
	}

	power /= 10;
	for (uint32 i = 0; power; ++i, power /= 10)
	{
		int32 d = value / power;
		str += (char)('0' + d);
		value -= d * power;
	}

	return str;
}


inline bool extend_string(uchar *&str, uint32 &size)
{
	uint32 new_size = size * 2;
	uchar *p = new uchar[new_size+1];

	if (!p)
		return false;

	std::copy(str, str+size, p);
	size = new_size;
	delete[] str;
	str = p;

	return true;
}


inline bool extend_string_to(uchar *&str, uint32 &size, uint32 new_size)
{
	if (new_size <= size)
		return true;

	uchar *p = new uchar[new_size+1];

	if (!p)
		return false;

	std::copy(str, str+size, p);

	size = new_size;
	delete[] str;
	str = p;

	return true;
}


inline uint32 int_log(uint64 x, uint64 base)
{
	uint32 r = 0;

	if (base == 0)
		return 1;
	if (base == 1)
		base++;

	for (uint64 tmp = base; tmp <= x; tmp *= base)
		++r;

	return r;
}


inline uint64 to_num(const uchar *str, uint32 len)
{
	uint64 r = 0;

	for (uint32 i = 0; i < len; ++i)
		r = r * 10 + (uint64)(str[i] - '0');

	return r;
}


inline bool is_num(const char* str_, uint32 len_, uint64& val_)
{
	val_ = 0;
	uint32 i;
	for (i = 0; i < len_; ++i)
	{
		if (str_[i] < '0' || str_[i] > '9')
			break;
		val_ = val_ * 10 + (str_[i] - '0');
	}

	return i == len_ && (len_ == 1 || str_[0] != '0');
}


inline uint32 bit_length(uint64 x)
{
	for (uint32 i = 0; i < 32; ++i)
	{
		if (x < (1ull << i))
			return i;
	}
	return 64;
}


#endif // H_UTILS
