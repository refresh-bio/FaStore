#include <stdio.h>
#include "rle.h"

int main(int argc_, char* argv_[])
{
	if (argc_ != 4 || (argv_[1][0] != 'e' && argv_[1][0] != 'd'))
	{
		fprintf(stderr, "usage: rle <e|d> <in_file> <out_file>\n");
		return RLE_ERR;
	}

	if (argv_[1][0] == 'e')
		return rle_encode(argv_[2], argv_[3]);
	return rle_decode(argv_[2], argv_[3]);
}
