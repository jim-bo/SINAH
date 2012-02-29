/*
 * types.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: eljimbo
 */

#include "types.h"

// node table
NodePair node_table_load(const char * file_path){

	// size.
	size_t NODETABLE_SIZE = sizeof( NodeTable );

	// offsets.
	size_t NODETABLE_OFFSET[NODE_TABLE_NFIELDS] = {
		HOFFSET( NodeTable, node_idx ),
		HOFFSET( NodeTable, ctg_name ),
		HOFFSET( NodeTable, ctg_width ),
		HOFFSET( NodeTable, ctg_orien ),
		HOFFSET( NodeTable, ctg_order ),
		HOFFSET( NodeTable, invalid )
	};

	// sizes.
	NodeTable NODETABLE_PROTO;
	size_t NODETABLE_SIZES[NODE_TABLE_NFIELDS] = {
		sizeof( NODETABLE_PROTO.node_idx ),
		sizeof( NODETABLE_PROTO.ctg_name ),
		sizeof( NODETABLE_PROTO.ctg_width ),
		sizeof( NODETABLE_PROTO.ctg_orien ),
		sizeof( NODETABLE_PROTO.ctg_order ),
		sizeof( NODETABLE_PROTO.invalid )
	};

	// open the file.
	hid_t file_id = H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);

	// get table description.
	hsize_t nfields;
	hsize_t nrecords;
	H5TBget_table_info(file_id, NODE_TABLE_NAME, &nfields, &nrecords );

	// read table into memory.
	NodeTable * buffer = new NodeTable[nrecords];
	H5TBread_table(file_id, NODE_TABLE_NAME, NODETABLE_SIZE, NODETABLE_OFFSET, NODETABLE_SIZES, buffer);

	// close out.
	H5Fclose( file_id );

	// return buffer array.
	return std::pair<NodeTable *, hsize_t>(buffer, nrecords);
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
	hid_t field_type[NODE_TABLE_NFIELDS];
	field_type[0] = H5T_NATIVE_LONG;		// node_idx
	field_type[1] = st_ctg_name;		// ctg_name[255]
	field_type[2] = H5T_NATIVE_LONG;		// ctg_width
	field_type[3] = H5T_NATIVE_LONG;		// ctg_orien
	field_type[4] = H5T_NATIVE_LONG;		// ctg_order
	field_type[5] = H5T_NATIVE_LONG;		// invalid

	// size.
	size_t NODETABLE_SIZE = sizeof( NodeTable );

	// offsets.
	size_t NODETABLE_OFFSET[NODE_TABLE_NFIELDS] = {
		HOFFSET( NodeTable, node_idx ),
		HOFFSET( NodeTable, ctg_name ),
		HOFFSET( NodeTable, ctg_width ),
		HOFFSET( NodeTable, ctg_orien ),
		HOFFSET( NodeTable, ctg_order ),
		HOFFSET( NodeTable, invalid )
	};

	// the rest.
	hid_t NODETABLE_FIELD_TYPE[NODE_TABLE_NFIELDS];
	hsize_t NODETABLE_CHUNK_SIZE = 10;
	const char * NODETABLE_FIELD_NAMES[NODE_TABLE_NFIELDS] = { "node_idx", "ctg_name", "ctg_width", "ctg_orien", "ctg_order", "invalid" };


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
		NODE_TABLE_NAME,
		NODE_TABLE_NFIELDS,
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

// edges
EdgePair edge_table_load(const char * file_path){

	// size.
	size_t EDGE_TABLE_SIZE = sizeof( EdgeTable );

	// offset.
	size_t EdgeTable_OFFSET[EDGE_TABLE_NFIELDS] = {
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

	// sizes.
	EdgeTable EdgeTable_PROTO;
	size_t EDGE_TABLE_SIZES[EDGE_TABLE_NFIELDS] = {
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

	// open the file.
	hid_t file_id = H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);

	// get table description.
	hsize_t nfields;
	hsize_t nrecords;
	H5TBget_table_info(file_id, EDGE_TABLE_NAME, &nfields, &nrecords );

	// read table into memory.
	EdgeTable * buffer = new EdgeTable[nrecords];
	H5TBread_table(file_id, EDGE_TABLE_NAME, EDGE_TABLE_SIZE, EdgeTable_OFFSET, EDGE_TABLE_SIZES, buffer);

	// close out.
	H5Fclose( file_id );

	// return buffer array.
	return EdgePair(buffer, nrecords);
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
	hid_t field_type[EDGE_TABLE_NFIELDS];
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

	// size.
	size_t EDGE_TABLE_SIZE = sizeof( EdgeTable );

	// offset.
	size_t EdgeTable_OFFSET[EDGE_TABLE_NFIELDS] = {
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

	// sizes.
	EdgeTable EdgeTable_PROTO;
	size_t EDGE_TABLE_SIZES[EDGE_TABLE_NFIELDS] = {
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

	// the rest.
	hid_t EdgeTable_FIELD_TYPE[EDGE_TABLE_NFIELDS];
	hsize_t EdgeTable_CHUNK_SIZE = 20;
	const char * EdgeTable_FIELD_NAMES[EDGE_TABLE_NFIELDS] = {"ctg_a_idx", "ctg_b_idx", "read_a_left_pos", "read_a_right_pos",
			"read_b_left_pos", "read_b_right_pos", "read_a_orien", "read_b_orien", "insert_size", "implied_state", "implied_dist", "std_dev", "invalid", "used"};

	// make the table.
	status = H5TBmake_table(
		table_title,
		file_id,
		EDGE_TABLE_NAME,
		EDGE_TABLE_NFIELDS,
		nrecords,
		EDGE_TABLE_SIZE,
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

// bundles.
BundlePair bundle_table_load(const char * file_path){

	// size.
	size_t BundleTable_SIZE = sizeof( BundleTable );

	// offset.
	size_t BundleTable_OFFSET[BUNDLE_TABLE_NFIELDS] = {
		HOFFSET( BundleTable, ctg_a_idx ),
		HOFFSET( BundleTable, ctg_b_idx ),
		HOFFSET( BundleTable, WT_A ),
		HOFFSET( BundleTable, WT_B ),
		HOFFSET( BundleTable, WT_C ),
		HOFFSET( BundleTable, WT_D )
	};

	// sizes.
	BundleTable BundleTable_PROTO;
	size_t BundleTable_SIZES[BUNDLE_TABLE_NFIELDS] = {
		sizeof( BundleTable_PROTO.ctg_a_idx ),
		sizeof( BundleTable_PROTO.ctg_b_idx ),
		sizeof( BundleTable_PROTO.WT_A ),
		sizeof( BundleTable_PROTO.WT_B ),
		sizeof( BundleTable_PROTO.WT_C ),
		sizeof( BundleTable_PROTO.WT_D )
	};

	// open the file.
	hid_t file_id = H5Fopen(file_path, H5F_ACC_RDWR, H5P_DEFAULT);

	// get table description.
	hsize_t nfields;
	hsize_t nrecords;
	H5TBget_table_info(file_id, BUNDLE_TABLE_NAME, &nfields, &nrecords );

	// read table into memory.
	BundleTable * buffer = new BundleTable[nrecords];
	H5TBread_table(file_id, BUNDLE_TABLE_NAME, BundleTable_SIZE, BundleTable_OFFSET, BundleTable_SIZES, buffer);

	// close out.
	H5Fclose( file_id );

	// return buffer array.
	return BundlePair(buffer, nrecords);
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
	hid_t field_type[BUNDLE_TABLE_NFIELDS];
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

	// size.
	size_t BundleTable_SIZE = sizeof( BundleTable );

	// offset.
	size_t BundleTable_OFFSET[BUNDLE_TABLE_NFIELDS] = {
		HOFFSET( BundleTable, ctg_a_idx ),
		HOFFSET( BundleTable, ctg_b_idx ),
		HOFFSET( BundleTable, WT_A ),
		HOFFSET( BundleTable, WT_B ),
		HOFFSET( BundleTable, WT_C ),
		HOFFSET( BundleTable, WT_D )
	};

	// sizes.
	BundleTable BundleTable_PROTO;
	size_t BundleTable_SIZES[BUNDLE_TABLE_NFIELDS] = {
		sizeof( BundleTable_PROTO.ctg_a_idx ),
		sizeof( BundleTable_PROTO.ctg_b_idx ),
		sizeof( BundleTable_PROTO.WT_A ),
		sizeof( BundleTable_PROTO.WT_B ),
		sizeof( BundleTable_PROTO.WT_C ),
		sizeof( BundleTable_PROTO.WT_D )
	};

	// the rest.
	hid_t BundleTable_FIELD_TYPE[BUNDLE_TABLE_NFIELDS];
	hsize_t BundleTable_CHUNK_SIZE = 10;
	const char * BundleTable_FIELD_NAMES[BUNDLE_TABLE_NFIELDS] = {"ctg_a_idx", "ctg_b_idx", "WT_A", "WT_B", "WT_C", "WT_D"};

	// make the table.
	status = H5TBmake_table(
		table_title,
		file_id,
		BUNDLE_TABLE_NAME,
		BUNDLE_TABLE_NFIELDS,
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



