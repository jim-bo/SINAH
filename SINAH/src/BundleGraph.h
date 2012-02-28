/*
 * BundleGraph.h
 *
 *  Created on: Feb 27, 2012
 *      Author: jlindsay
 */

// include OGDF
#include "OGDF.h"

#ifndef BUNDLEGRAPH_H_
#define BUNDLEGRAPH_H_

class BundleGraph {
public:
	// the graph.
	ogdf::Graph G;

	// index to node/edge vectors.
	vector<ogdf::node> idx2n;

	// node/edge to idnex lookup array.
	ogdf::NodeArray<int> n2idx;

}


#endif /* BUNDLEGRAPH_H_ */
