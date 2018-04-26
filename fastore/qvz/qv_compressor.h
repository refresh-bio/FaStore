/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on QVZ software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef qv_compressor_h
#define qv_compressor_h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codebook.h"
#include "../core/BitMemory.h"

#define m_arith  22

#define OS_STREAM_BUF_LEN		(4096*4096)

#define COMPRESSION 0
#define DECOMPRESSION 1

typedef struct Arithmetic_code_t {
    int32_t scale3;
    
	uint32_t l;
    uint32_t u;
    uint32_t t;

    uint32_t m;
	uint32_t r;			// Rescaling condition
}*Arithmetic_code;


typedef struct os_stream_t {
	//FILE *fp;
	uint8_t *buf;
	uint32_t bufPos;
	uint8_t bitPos;
	uint64_t written;
} *osStream;

typedef struct stream_stats_t {
    uint32_t *counts;
    uint32_t alphabetCard;
    uint32_t step;
    uint32_t n;
} *stream_stats_ptr_t;

typedef struct arithStream_t {
	stream_stats_ptr_t cluster_stats;
    stream_stats_ptr_t ***stats;
	Arithmetic_code a;
	BitMemoryReader* is;
	BitMemoryWriter* os;
}*arithStream;

typedef struct qv_compressor_t{
    arithStream Quals;
}*qv_compressor;




// Stream interface
/*
struct os_stream_t *alloc_os_stream(FILE *fp, uint8_t in);
void free_os_stream(struct os_stream_t *);
uint8_t stream_read_bit(struct os_stream_t *);
uint32_t stream_read_bits(struct os_stream_t *os, uint8_t len);
void stream_write_bit(struct os_stream_t *, uint8_t);
void stream_write_bits(struct os_stream_t *os, uint32_t dw, uint8_t len);
void stream_finish_byte(struct os_stream_t *);
void stream_write_buffer(struct os_stream_t *);
*/

inline uint8_t stream_read_bit(BitMemoryReader* is)
{
	return (uint8_t)is->GetBit();
}

inline uint32_t stream_read_bits(BitMemoryReader *is, uint8_t len)
{
	return is->GetBits(len);
}

inline void stream_write_bit(BitMemoryWriter *os, uint8_t b)
{
	os->PutBit(b);
}

inline void stream_write_bits(BitMemoryWriter *os, uint32_t dw, uint8_t len)
{
	os->PutBits(dw, len);
}

inline void stream_finish_byte(struct BitMemoryWriter *os)
{
	os->FillLastByte();
}

// Arithmetic ncoder interface
Arithmetic_code initialize_arithmetic_encoder(uint32_t m);
void arithmetic_encoder_step(Arithmetic_code a, stream_stats_ptr_t stats, int32_t x, BitMemoryWriter* os);
uint64_t encoder_last_step(Arithmetic_code a, BitMemoryWriter* os);
uint32_t arithmetic_decoder_step(Arithmetic_code a, stream_stats_ptr_t stats, BitMemoryReader* is);
uint32_t decoder_last_step(Arithmetic_code a, stream_stats_ptr_t stats);

// Encoding stats management
stream_stats_ptr_t **initialize_stream_stats(struct cond_quantizer_list_t *q_list);
void update_stats(stream_stats_ptr_t stats, uint32_t x, uint32_t r);

// Quality value compression interface
void compress_qv(arithStream as, uint32_t x, uint8_t cluster, uint32_t column, uint32_t idx);
void qv_write_cluster(arithStream as, uint8_t cluster);
uint32_t decompress_qv(arithStream as, uint8_t cluster, uint32_t column, uint32_t idx);
uint8_t qv_read_cluster(arithStream as);

bool initialize_qvz_compressor(arithStream as, cond_quantizer_list_t *q);
bool initialize_qvz_decompressor(arithStream as, cond_quantizer_list_t *q);
void free_qvz_arith(arithStream as, struct cond_quantizer_list_t *info);

uint64_t start_qv_compression(struct quality_file_t *info, FILE *fout, double *dis, FILE * funcompressed);
void start_qv_decompression(FILE *fout, FILE *fin, struct quality_file_t *info);


// QVZ encoder <---> ORCOM encoder : proxy
//
class QVZCoder
{
public:
	QVZCoder(cond_quantizer_list_t *q)
		:	q_info(q)
	{
		qv = (qv_compressor)calloc(1, sizeof(struct qv_compressor_t));
		as = (arithStream)calloc(1, sizeof(struct arithStream_t));
		qv->Quals = as;
	}

	virtual ~QVZCoder()
	{
		free(qv);
		free(as);
	}

protected:
	qv_compressor qv;
	arithStream as;
	cond_quantizer_list_t *q_info;
};

class QVZEncoder : public QVZCoder
{
public:
	QVZEncoder(BitMemoryWriter* writer, cond_quantizer_list_t *q)
		:	QVZCoder(q)
	{
		as->os = writer;

		initialize_qvz_compressor(as, q);
	}

	~QVZEncoder()
	{
		free_qvz_arith(as, q_info);
	}


	void Start()
	{
		as->a = initialize_arithmetic_encoder(m_arith);
		as->a->t = 0;
	}

	void EncodeNext(uint32 q_state, uint32 col, uint32 idx)
	{
		compress_qv(qv->Quals, q_state, 0, col, idx);
	}

	void End()
	{
		encoder_last_step(qv->Quals->a, qv->Quals->os);

		free(as->a);
	}

private:

};

class QVZDecoder : public QVZCoder
{
public:
	QVZDecoder(BitMemoryReader* reader, cond_quantizer_list_t *q)
		:	QVZCoder(q)
	{
		as->is = reader;

		initialize_qvz_decompressor(as, q);
	}

	~QVZDecoder()
	{
		free_qvz_arith(as, q_info);
	}


	void Start()
	{
		as->a = initialize_arithmetic_encoder(m_arith);
		as->a->t = as->is->GetBits(as->a->m);
	}

	uint32 DecodeNext(uint32 col, uint32 idx)
	{
		return decompress_qv(qv->Quals, 0, col, idx);
	}

	void End()
	{
		free(as->a);
	}

private:
};

#endif
