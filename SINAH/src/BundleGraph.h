/*
 * BundleGraph.h
 *
 *  Created on: Feb 27, 2012
 *      Author: jlindsay
 */

#ifndef BUNDLEGRAPH_H_
#define BUNDLEGRAPH_H_

// google
#include <google/sparse_hash_map>
#include <google/sparse_hash_set>

// include OGDF
#include "OGDF.h"

// include types.
#include "types.h"

// define a node set.
typedef google::sparse_hash_set<int> NodeSet;

class BundleGraph {
public:
	// the graph.
	ogdf::Graph G;

	// index to node map.
	google::sparse_hash_map<int,ogdf::node> idx2n;

	// node/edge to index lookup array.
	ogdf::NodeArray<int> n2idx;

	// constructor.
	BundleGraph(NodePair np, BundlePair bp, NodeSet active);

};


#endif /* BUNDLEGRAPH_H_ */
