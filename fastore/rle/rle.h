#ifndef H_RLE
#define H_RLE

#include <stdint.h>

#define RLE_OK 0
#define RLE_ERR -1

int rle_encode(const char* in_file_, const char* out_file_);
int rle_encode(unsigned char* in_mem_, uint64_t in_size_, unsigned char* out_mem_, uint64_t* out_size_);

int rle_decode(const char* in_file_, const char* out_file_);
int rle_decode(unsigned char* in_mem_, uint64_t in_size_, unsigned char* out_mem_, uint64_t* out_size_);


#endif // H_RLE
