#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

// necessary libraries.
#include "hdf5.h"
#include "hdf5_hl.h"
using namespace std;

#define AGP_TABLE_NAME "AgpTable"
#define AGP_NFIELDS  (hsize_t)  10
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
const char * AGP_FIELD_NAMES[AGP_NFIELDS]  = {
	"scaf_name",
	"scaf_start", 
	"scaf_stop", 
	"scaf_idx", 
	"comp_type", 
	"comp_name", 
	"comp_start", 
	"comp_stop", 
	"comp_orien", 
	"comp_linkage"
};
int save_agp(int num_records, char * file_name, char * table_title, AgpTable * p_data){
	
	// setup table variables.
	size_t NRECORDS = (size_t) num_records;
	size_t AGP_SIZE = sizeof(AgpTable);

	size_t AGP_OFFSET[AGP_NFIELDS] = { 
		HOFFSET( AgpTable, scaf_name ),
		HOFFSET( AgpTable, scaf_start ),
		HOFFSET( AgpTable, scaf_stop ),
		HOFFSET( AgpTable, scaf_idx ),
		HOFFSET( AgpTable, comp_type ),
		HOFFSET( AgpTable, comp_name ),
		HOFFSET( AgpTable, comp_start ),
		HOFFSET( AgpTable, comp_stop ),
		HOFFSET( AgpTable, comp_orien ),
		HOFFSET( AgpTable, comp_linkage )
	};

	size_t AGP_SIZES[AGP_NFIELDS] = { 
		sizeof( p_data[0].scaf_name ),
		sizeof( p_data[0].scaf_start ),
		sizeof( p_data[0].scaf_stop ),
		sizeof( p_data[0].scaf_idx ),
		sizeof( p_data[0].comp_type ),
		sizeof( p_data[0].comp_name ),
		sizeof( p_data[0].comp_start ),
		sizeof( p_data[0].comp_stop ),
		sizeof( p_data[0].comp_orien ),
		sizeof( p_data[0].comp_linkage )
	};
	
	hid_t field_type[AGP_NFIELDS];
	hid_t file_id;
	hsize_t chunk_size = 10;
	herr_t status;
	int * fill_data = NULL;
	int compress  = 0;
	int i;
	
	// open string types.
	hid_t st_scaf_name = H5Tcopy( H5T_C_S1 );
	hid_t st_comp_type = H5Tcopy( H5T_C_S1 );
	hid_t st_comp_name = H5Tcopy( H5T_C_S1 );
	
	H5Tset_size( st_scaf_name, 255 );
	H5Tset_size( st_comp_type, 255 );
	H5Tset_size( st_comp_name, 255 );
	
	// save field types.
	field_type[0] = st_scaf_name;		// scaf_name
	field_type[1] = H5T_NATIVE_LONG;		// scaf_start
	field_type[2] = H5T_NATIVE_LONG;		// scaf_stop
	field_type[3] = H5T_NATIVE_LONG;		// scaf_idx
	field_type[4] = st_comp_type;		// comp_type
	field_type[5] = st_comp_name;		// comp_name
	field_type[6] = H5T_NATIVE_LONG;		// comp_start
	field_type[7] = H5T_NATIVE_LONG;		// comp_stop
	field_type[8] = H5T_NATIVE_LONG;		// comp_orien
	field_type[9] = H5T_NATIVE_LONG;		// comp_linkage
	
	// create new file.
	file_id = H5Fcreate( 
		file_name, 
		H5F_ACC_TRUNC, 
		H5P_DEFAULT, 
		H5P_DEFAULT 
	);

	// make the table.
	status = H5TBmake_table( 
		table_title, 
		file_id, 
		AGP_TABLE_NAME,
		AGP_NFIELDS,
		NRECORDS,
		AGP_SIZE,
		AGP_FIELD_NAMES,
		AGP_OFFSET,
		field_type,
		chunk_size, 
		fill_data, 
		compress, 
		p_data  
	);

	// close stirngs.
	H5Tclose( st_scaf_name );
	H5Tclose( st_comp_type );
	H5Tclose( st_comp_name );
 
 /* close the file */
 H5Fclose( file_id );


	// return file_id.
	return file_id;	
}

//------------------ NODE TABLE -----------------//
//------------------ NODE TABLE -----------------//
//------------------ NODE TABLE -----------------//
typedef struct NodeTable {
	long    node_idx;
	char   	ctg_name[255];
	long   	ctg_width;
	long   	ctg_orien;
	long   	ctg_order;
	long	invalid;
} NodeTable;
NodeTable NODETABLE_PROTO;
#define NODETABLE_NFIELDS  (hsize_t)  6
#define NODETABLE_NAME "nodes"
size_t NODETABLE_SIZE = sizeof( NodeTable );
size_t NODETABLE_OFFSET[NODETABLE_NFIELDS] = {
	HOFFSET( NodeTable, node_idx ),
	HOFFSET( NodeTable, ctg_name ),
	HOFFSET( NodeTable, ctg_width ),
	HOFFSET( NodeTable, ctg_orien ),
	HOFFSET( NodeTable, ctg_order ),
	HOFFSET( NodeTable, invalid )
};
size_t NODETABLE_SIZES[NODETABLE_NFIELDS] = {
	sizeof( NODETABLE_PROTO.node_idx ),
	sizeof( NODETABLE_PROTO.ctg_name ),
	sizeof( NODETABLE_PROTO.ctg_width ),
	sizeof( NODETABLE_PROTO.ctg_orien ),
	sizeof( NODETABLE_PROTO.ctg_order ),
	sizeof( NODETABLE_PROTO.invalid )
};
hid_t NODETABLE_FIELD_TYPE[NODETABLE_NFIELDS];
hsize_t NODETABLE_CHUNK_SIZE = 10;
const char * NODETABLE_FIELD_NAMES[NODETABLE_NFIELDS] = { "node_idx", "ctg_name", "ctg_width", "ctg_orien", "ctg_order", "invalid" };

pair<NodeTable *, hsize_t> node_table_load(const char * file_path){

	// open the file.
	hid_t file_id = H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);
	
	// get table description.
	hsize_t nfields;
	hsize_t nrecords;
	H5TBget_table_info(file_id, NODETABLE_NAME, &nfields, &nrecords );
	
	// read table into memory.
	NodeTable * buffer = new NodeTable[nrecords];
	H5TBread_table(file_id, NODETABLE_NAME, NODETABLE_SIZE, NODETABLE_OFFSET, NODETABLE_SIZES, buffer);

	// close out.
	H5Fclose( file_id );

	// return buffer array.
	return pair<NodeTable *, hsize_t>(buffer, nrecords);
}
void node_table_save(NodeTable * buffer, int nr, const char * table_title, const char * file_path){

	// create basic info.
	hsize_t nrecords = (hsize_t) nr;
	hid_t file_id;
	herr_t status;
	int * fill_data = NULL;
	int compress  = 0;
	int i;
	
	// open string types.
	hid_t st_ctg_name = H5Tcopy( H5T_C_S1 );
	H5Tset_size( st_ctg_name, 255 );
	
	// save field types.
	hid_t field_type[NODETABLE_NFIELDS];
	field_type[0] = H5T_NATIVE_LONG;		// node_idx
	field_type[1] = st_ctg_name;		// ctg_name[255]
	field_type[2] = H5T_NATIVE_LONG;		// ctg_width
	field_type[3] = H5T_NATIVE_LONG;		// ctg_orien
	field_type[4] = H5T_NATIVE_LONG;		// ctg_order
	field_type[5] = H5T_NATIVE_LONG;		// invalid
	
	// create new file.
	file_id = H5Fcreate( 
		file_path, 
		H5F_ACC_TRUNC, 
		H5P_DEFAULT, 
		H5P_DEFAULT 
	);

	// make the table.
	status = H5TBmake_table( 
		table_title, 
		file_id, 
		NODETABLE_NAME,
		NODETABLE_NFIELDS,
		nrecords,
		NODETABLE_SIZE,
		NODETABLE_FIELD_NAMES,
		NODETABLE_OFFSET,
		field_type,
		NODETABLE_CHUNK_SIZE, 
		fill_data, 
		compress, 
		buffer  
	);

	// close strings.
	H5Tclose( st_ctg_name );
 
	// close name.
	H5Fclose( file_id );
}

//------------------ EDGE TABLE -----------------//
//------------------ EDGE TABLE -----------------//
//------------------ EDGE TABLE -----------------//
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
EdgeTable EdgeTable_PROTO;
#define EdgeTable_NFIELDS  (hsize_t)  14
#define EdgeTable_NAME "edges"
size_t EdgeTable_SIZE = sizeof( EdgeTable );
size_t EdgeTable_OFFSET[EdgeTable_NFIELDS] = {
	HOFFSET( EdgeTable, ctg_a_idx ),
	HOFFSET( EdgeTable, ctg_b_idx ),
	HOFFSET( EdgeTable, read_a_left_pos ),
	HOFFSET( EdgeTable, read_a_right_pos ),
	HOFFSET( EdgeTable, read_b_left_pos ),
	HOFFSET( EdgeTable, read_b_right_pos ),
	HOFFSET( EdgeTable, read_a_orien ),
	HOFFSET( EdgeTable, read_b_orien ),
	HOFFSET( EdgeTable, insert_size ),
	HOFFSET( EdgeTable, implied_state ),
	HOFFSET( EdgeTable, implied_dist ),
	HOFFSET( EdgeTable, std_dev ),
	HOFFSET( EdgeTable, invalid ),
	HOFFSET( EdgeTable, used )
};
size_t EdgeTable_SIZES[EdgeTable_NFIELDS] = {
	sizeof( EdgeTable_PROTO.ctg_a_idx ),
	sizeof( EdgeTable_PROTO.ctg_b_idx ),
	sizeof( EdgeTable_PROTO.read_a_left_pos ),
	sizeof( EdgeTable_PROTO.read_a_right_pos ),
	sizeof( EdgeTable_PROTO.read_b_left_pos ),
	sizeof( EdgeTable_PROTO.read_b_right_pos ),	
	sizeof( EdgeTable_PROTO.read_a_orien ),	
	sizeof( EdgeTable_PROTO.read_b_orien ),	
	sizeof( EdgeTable_PROTO.insert_size ),	
	sizeof( EdgeTable_PROTO.implied_state ),	
	sizeof( EdgeTable_PROTO.implied_dist ),	
	sizeof( EdgeTable_PROTO.std_dev ),	
	sizeof( EdgeTable_PROTO.invalid ),
	sizeof( EdgeTable_PROTO.used )
};
hid_t EdgeTable_FIELD_TYPE[EdgeTable_NFIELDS];
hsize_t EdgeTable_CHUNK_SIZE = 20;
const char * EdgeTable_FIELD_NAMES[EdgeTable_NFIELDS] = {"ctg_a_idx", "ctg_b_idx", "read_a_left_pos", "read_a_right_pos",
	"read_b_left_pos", "read_b_right_pos", "read_a_orien", "read_b_orien", "insert_size", "implied_state", "implied_dist", "std_dev", "invalid", "used"};

pair<EdgeTable *, hsize_t> edge_table_load(const char * file_path){

	// open the file.
	hid_t file_id = H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);
	
	// get table description.
	hsize_t nfields;
	hsize_t nrecords;
	H5TBget_table_info(file_id, EdgeTable_NAME, &nfields, &nrecords );
	
	// read table into memory.
	EdgeTable * buffer = new EdgeTable[nrecords];
	H5TBread_table(file_id, EdgeTable_NAME, EdgeTable_SIZE, EdgeTable_OFFSET, EdgeTable_SIZES, buffer);

	// close out.
	H5Fclose( file_id );

	// return buffer array.
	return pair<EdgeTable *, hsize_t>(buffer, nrecords);
}
void edge_table_save(EdgeTable * buffer, int nr, const char * table_title, const char * file_path){

	// create basic info.
	hsize_t nrecords = (hsize_t) nr;
	hid_t file_id;
	herr_t status;
	int * fill_data = NULL;
	int compress  = 0;
	int i;
	
	// save field types.
	hid_t field_type[EdgeTable_NFIELDS];
	field_type[0] = H5T_NATIVE_LONG;		// node_idx
	field_type[1] = H5T_NATIVE_LONG;		// node_idx
	field_type[2] = H5T_NATIVE_LONG;		// ctg_width
	field_type[3] = H5T_NATIVE_LONG;		// ctg_orien
	field_type[4] = H5T_NATIVE_LONG;		// ctg_order
	field_type[5] = H5T_NATIVE_LONG;		// invalid
	field_type[6] = H5T_NATIVE_LONG;		// invalid
	field_type[7] = H5T_NATIVE_LONG;		// invalid
	field_type[8] = H5T_NATIVE_LONG;		// invalid
	field_type[9] = H5T_NATIVE_LONG;		// invalid
	field_type[10] = H5T_NATIVE_LONG;	// invalid
	field_type[11] = H5T_NATIVE_LONG;	// invalid
	field_type[12] = H5T_NATIVE_LONG;	// invalid
	field_type[13] = H5T_NATIVE_LONG;	// invalid
	
	// create new file.
	file_id = H5Fcreate( 
		file_path, 
		H5F_ACC_TRUNC, 
		H5P_DEFAULT, 
		H5P_DEFAULT 
	);

	// make the table.
	status = H5TBmake_table( 
		table_title, 
		file_id, 
		EdgeTable_NAME,
		EdgeTable_NFIELDS,
		nrecords,
		EdgeTable_SIZE,
		EdgeTable_FIELD_NAMES,
		EdgeTable_OFFSET,
		field_type,
		EdgeTable_CHUNK_SIZE, 
		fill_data, 
		compress, 
		buffer  
	);

	// close name.
	H5Fclose( file_id );
}

//------------------ BUNDLE TABLE -----------------//
//------------------ BUNDLE TABLE -----------------//
//------------------ BUNDLE TABLE -----------------//
typedef struct BundleTable {
	long    	ctg_a_idx;
	long    	ctg_b_idx;
	double  	WT_A;
	double   	WT_B;
	double   	WT_C;
	double   	WT_D;
} BundleTable;
BundleTable BundleTable_PROTO;
#define BundleTable_NFIELDS  (hsize_t)  6
#define BundleTable_NAME "bundles"
size_t BundleTable_SIZE = sizeof( BundleTable );
size_t BundleTable_OFFSET[BundleTable_NFIELDS] = {
	HOFFSET( BundleTable, ctg_a_idx ),
	HOFFSET( BundleTable, ctg_b_idx ),
	HOFFSET( BundleTable, WT_A ),
	HOFFSET( BundleTable, WT_B ),
	HOFFSET( BundleTable, WT_C ),
	HOFFSET( BundleTable, WT_D )
};
size_t BundleTable_SIZES[BundleTable_NFIELDS] = {
	sizeof( BundleTable_PROTO.ctg_a_idx ),
	sizeof( BundleTable_PROTO.ctg_b_idx ),
	sizeof( BundleTable_PROTO.WT_A ),
	sizeof( BundleTable_PROTO.WT_B ),
	sizeof( BundleTable_PROTO.WT_C ),
	sizeof( BundleTable_PROTO.WT_D )		
};
hid_t BundleTable_FIELD_TYPE[BundleTable_NFIELDS];
hsize_t BundleTable_CHUNK_SIZE = 10;
const char * BundleTable_FIELD_NAMES[BundleTable_NFIELDS] = {"ctg_a_idx", "ctg_b_idx", "WT_A", "WT_B", 	"WT_C", "WT_D"};

pair<BundleTable *, hsize_t> bundle_table_load(const char * file_path){

	// open the file.
	hid_t file_id = H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);
	
	// get table description.
	hsize_t nfields;
	hsize_t nrecords;
	H5TBget_table_info(file_id, BundleTable_NAME, &nfields, &nrecords );
	
	// read table into memory.
	BundleTable * buffer = new BundleTable[nrecords];
	H5TBread_table(file_id, BundleTable_NAME, BundleTable_SIZE, BundleTable_OFFSET, BundleTable_SIZES, buffer);

	// close out.
	H5Fclose( file_id );

	// return buffer array.
	return pair<BundleTable *, hsize_t>(buffer, nrecords);
}
void bundle_table_save(BundleTable * buffer, int nr, const char * table_title, const char * file_path){

	// create basic info.
	hsize_t nrecords = (hsize_t) nr;
	hid_t file_id;
	herr_t status;
	int * fill_data = NULL;
	int compress  = 0;
	int i;
	
	// save field types.
	hid_t field_type[BundleTable_NFIELDS];
	field_type[0] = H5T_NATIVE_LONG;		// node_idx
	field_type[1] = H5T_NATIVE_LONG;		// node_idx
	field_type[2] = H5T_NATIVE_DOUBLE;		// ctg_width
	field_type[3] = H5T_NATIVE_DOUBLE;		// ctg_orien
	field_type[4] = H5T_NATIVE_DOUBLE;		// ctg_order
	field_type[5] = H5T_NATIVE_DOUBLE;		// invalid
	
	// create new file.
	file_id = H5Fcreate( 
		file_path, 
		H5F_ACC_TRUNC, 
		H5P_DEFAULT, 
		H5P_DEFAULT 
	);

	// make the table.
	status = H5TBmake_table( 
		table_title, 
		file_id, 
		BundleTable_NAME,
		BundleTable_NFIELDS,
		nrecords,
		BundleTable_SIZE,
		BundleTable_FIELD_NAMES,
		BundleTable_OFFSET,
		field_type,
		BundleTable_CHUNK_SIZE, 
		fill_data, 
		compress, 
		buffer  
	);

	// close name.
	H5Fclose( file_id );
}

/* typedef section */
typedef pair<AgpTable *, hsize_t> AgpPair;
typedef pair<NodeTable *, hsize_t> NodePair;
typedef pair<EdgeTable *, hsize_t> EdgePair;
typedef pair<BundleTable *, hsize_t> BundlePair;

#endif
