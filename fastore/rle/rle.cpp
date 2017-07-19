#include <stdio.h>
#include <stdint.h>

#include "rle.h"

#define RLE_OFFSET 2


/***************
*
* file-to-file encode
*
*/
int rle_encode(const char* in_file_, const char* out_file_)
{
	FILE* in = fopen(in_file_, "rb");
	if (in == NULL)
	{
		fprintf(stderr, "Cannot open input file: %s\n", in_file_);
		return RLE_ERR;
	}

	FILE* out = fopen(out_file_, "wb");
	if (out == NULL)
	{
		fprintf(stderr, "Cannot open output file: %s\n", out_file_);
		fclose(in);
		return RLE_ERR;
	}

	int cnt = 0;
	int c;
	while ((c = getc(in)) != EOF)
	{
		if (c == '.')
		{
			cnt++;
			if(cnt == 255 - RLE_OFFSET)
			{
				putc(cnt + RLE_OFFSET, out);
				cnt = 0;
			}
		}
		else
		{
			bool must_be_mis = (cnt > 0) && (cnt < 255 - RLE_OFFSET);
			if (cnt > 0)
			{
				putc(cnt + RLE_OFFSET, out);
				cnt = 0;
			}
			if (!must_be_mis)
				putc(0, out);
		}	
	}
	if (cnt > 0)
		putc(cnt + RLE_OFFSET, out);

	fclose(in);
	fclose(out);
	return RLE_OK;
}


/***************
*
* in-memory encode
*
*/
int rle_encode(unsigned char* in_mem_, uint64_t in_size_, unsigned char* out_mem_, uint64_t* out_size_)
{
	int cnt = 0;
	int64_t s = 0;
	*out_size_ = 0;
	for (uint64_t i = 0; i < in_size_; ++i)
	{
		int c = in_mem_[i];
		if (c == '.')
		{
			cnt++;
			if(cnt == 255 - RLE_OFFSET)
			{
				out_mem_[s++] = cnt + RLE_OFFSET;
				cnt = 0;
			}
		}
		else
		{
			bool must_be_mis = (cnt > 0) && (cnt < 255 - RLE_OFFSET);
			if (cnt > 0)
			{
				out_mem_[s++] = cnt + RLE_OFFSET;
				cnt = 0;
			}
			if (!must_be_mis)
				out_mem_[s++] = 0;
		}	
	}
	if (cnt > 0)
		out_mem_[s++] = cnt + RLE_OFFSET;

	*out_size_ = s;
	return RLE_OK;
}


/***************
*
* file-to-file decode
*
*/
int rle_decode(const char* in_file_, const char* out_file_)
{
	FILE* in = fopen(in_file_, "rb");
	if (in == NULL)
	{
		fprintf(stderr, "Cannot open input file: %s\n", in_file_);
		return RLE_ERR;
	}

	FILE* out = fopen(out_file_, "wb");
	if (out == NULL)
	{
		fprintf(stderr, "Cannot open output file: %s\n", out_file_);
		fclose(in);
		return RLE_ERR;
	}
	
	int c;
	while ((c = getc(in)) != EOF)
	{
		if (c > 0)
		{
			for (int i = 0; i < c - RLE_OFFSET; ++i)
				putc('.', out);

			if (c < 255)
				putc('*', out);
		}
		else
			putc('*', out);
	}
			
	fclose(in);
	fclose(out);
	return RLE_OK;
}


/***************
*
* in-memory decode
*
*/
int rle_decode(unsigned char* in_mem_, uint64_t in_size_, unsigned char* out_mem_, uint64_t* out_size_)
{
	uint64_t s = 0;
	for (uint64_t i = 0; i < in_size_; ++i)
	{
		int cnt = in_mem_[i];
		if (cnt > 0)
		{
			for (int j = 0; j < cnt - RLE_OFFSET; ++j)
				out_mem_[s++] = '.';

			if (cnt < 255)
				out_mem_[s++] = '*';
		}
		else
			out_mem_[s++] = '*';
	}
	*out_size_ = s;
	return RLE_OK;
}
