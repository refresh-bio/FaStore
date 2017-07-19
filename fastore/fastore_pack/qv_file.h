#ifndef _LINES_H_
#define _LINES_H_


#include <stdint.h>

#include "util.h"

#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <time.h>


#include "pmf.h"
#include "distortion.h"
#include "well.h"

#define ALPHABET_SIZE 72

#define qv2ch(a) a - 33

// This limits us to chunks that aren't too big to fit into a modest amount of memory at a time
#define MAX_LINES_PER_BLOCK			1000000
#define MAX_READS_PER_LINE			1022
#define READ_LINEBUF_LENGTH			(MAX_READS_PER_LINE+2)

// Error codes for reading a line block
#define LF_ERROR_NONE				0
#define LF_ERROR_NOT_FOUND			1
#define LF_ERROR_NO_MEMORY			2
#define LF_ERROR_TOO_LONG			4

/**
 * Points to a single line, which may be a pointer to a file in memory
 */
struct line_t {
	uint8_t cluster;		// Assigned cluster ID
	const symbol_t *m_data;	// Pointer to part of mmap'd region, has no offsets applied, do not modify!
};

/**
 * Points to a block of lines for incremental processing
 */
struct line_block_t {
	uint32_t count;
	struct line_t *lines;
};

/**
 * Stores information about a specific cluster
 */
struct cluster_t {
	// Used to do clustering
	uint8_t id;					// Cluster ID
	uint32_t count;				// Number of lines in this cluster
	//symbol_t *mean;				// Mean values for this cluster
	//uint64_t *accumulator;		// Accumulator for finding a new cluster center

	// Used after clustering is done
	struct cond_pmf_list_t *training_stats;
	struct cond_quantizer_list_t *qlist;
};

/**
 * Stores all clusters
 */
struct cluster_list_t {
	uint8_t count;
	struct cluster_t *clusters;
	//double *distances;			// Temporary storage for distances to each cluster center
};

/**
 *
 */
typedef struct qv_file_t{
    uint8_t* file_head;
    uint64_t lines;
    uint32_t read_length;
    uint32_t alphabet_size;
} *qv_file;

/**
 * Points to a file descriptor that includes important metadata about the file
 */
struct quality_file_t {
	struct alphabet_t *alphabet;
	//char *path;
	//uint64_t lines;
	uint32_t columns;
	//uint32_t block_count;
	//struct line_block_t *blocks;
	//uint8_t cluster_count;
	//struct cluster_list_t *clusters;
    struct cond_pmf_list_t *training_stats;
	struct distortion_t *dist;
	struct qv_options_t *opts;
    struct well_state_t well;
    struct cond_quantizer_list_t *qlist;
    //struct qv_file_t *qv_f;
};



// Memory management
uint32_t load_qv_file(const char *path, struct quality_file_t *info, uint64_t max_lines);
uint32_t alloc_blocks(struct quality_file_t *info);
void free_blocks(struct quality_file_t *info);

uint32_t generate_qv_struct(struct qv_file_t* qvs, struct quality_file_t *info, uint64_t max_lines);
qv_file load_file(const char *path, uint64_t max_lines);

#endif
