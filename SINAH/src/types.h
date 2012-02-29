#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

// necessary libraries.
#include "hdf5.h"
#include "hdf5_hl.h"
#include <utility>

// table names.
#define NODE_TABLE_NAME "nodes"
#define EDGE_TABLE_NAME "edges"
#define BUNDLE_TABLE_NAME "bundles"
#define AGP_TABLE_NAME "agps"

// table sizes.
#define NODE_TABLE_NFIELDS  (hsize_t)  6
#define EDGE_TABLE_NFIELDS  (hsize_t)  14
#define BUNDLE_TABLE_NFIELDS  (hsize_t)  6
#define AGP_TABLE_NFIELDS  (hsize_t)  10

// structs.
typedef struct NodeTable {
	long    node_idx;
	char   	ctg_name[255];
	long   	ctg_width;
	long   	ctg_orien;
	long   	ctg_order;
	long	invalid;
} NodeTable;

typedef struct EdgeTable {
	long    	ctg_a_idx;
	long    	ctg_b_idx;
	long    	read_a_left_pos;
	long    	read_a_right_pos;
	long    	read_b_left_pos;
	long    	read_b_right_pos;
	long    	read_a_orien;
	long    	read_b_orien;
	long    	insert_size;
	long    	implied_state;
	long    	implied_dist;
	long    	std_dev;
	long   		invalid;
	long   		used;
} EdgeTable;

typedef struct BundleTable {
	long    	ctg_a_idx;
	long    	ctg_b_idx;
	double  	WT_A;
	double   	WT_B;
	double   	WT_C;
	double   	WT_D;
} BundleTable;

typedef struct AgpTable {
	char   	scaf_name[255];
	long    	scaf_start;
	long		scaf_stop;
	long    	scaf_idx;
	char   	comp_type[255];
	char   	comp_name[255];
	long   	comp_start;
	long   	comp_stop;
	long   	comp_orien;
	long   	comp_linkage;
} AgpTable;

// pairs.
typedef std::pair<NodeTable *, hsize_t> NodePair;
typedef std::pair<EdgeTable *, hsize_t> EdgePair;
typedef std::pair<BundleTable *, hsize_t> BundlePair;
typedef std::pair<AgpTable *, hsize_t> AgpPair;

// load functions.
NodePair node_table_load(const char * file_path);
EdgePair edge_table_load(const char * file_path);
BundlePair bundle_table_load(const char * file_path);

// save functions.
void node_table_save(NodeTable * buffer, int nr, const char * table_title, const char * file_path);
void edge_table_save(EdgeTable * buffer, int nr, const char * table_title, const char * file_path);
void bundle_table_save(BundleTable * buffer, int nr, const char * table_title, const char * file_path);

#endif
