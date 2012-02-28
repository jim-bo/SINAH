/*
 * BundleGraph.h
 *
 *  Created on: Feb 27, 2012
 *      Author: jlindsay
 */

#ifndef BUNDLEGRAPH_H_
#define BUNDLEGRAPH_H_

// include master.
#include "SINAH.h"

// include OGDF
#include "OGDF.h"

// include types
#include "types.h"

// google
#include <google/sparse_hash_map>
#include <google/sparse_hash_set>

// define a node set.
typedef google::sparse_hash_set<int> NodeSet;

class BundleGraph {
public:
	// the graph.
	ogdf::Graph G;

	// index to node map.
	google::sparse_hash_map<int,ogdf::node> idx2n;

	// node/edge to idnex lookup array.
	ogdf::NodeArray<int> n2idx;

	// constructor.
	BundleGraph(NodePair np, BundlePair bp, NodeSet active);

};


#endif /* BUNDLEGRAPH_H_ */
