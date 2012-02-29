/*
 * BundleGraph.cpp
 *
 *  Created on: Feb 28, 2012
 *      Author: jlindsay
 */

// includes
#include "BundleGraph.h"

// namepsace
using namespace ogdf;

// basic constructor.
BundleGraph::BundleGraph(NodePair np, BundlePair bp, NodeSet active){

	// associate Node/Edge Array with graph.
	n2idx.init(G);

	// clear node map.
	idx2n.clear();

	// add all nodes.
	for(hsize_t i=0; i<np.second; i++){

		// skip check.
		if( active.find(np.first[i].node_idx) == active.end() ) continue;

		// add node.
		ogdf::node n = G.newNode();

		// track it.
		idx2n[np.first[i].node_idx] = n;
		n2idx[n] = np.first[i].node_idx;
	}

	// add all appropriate bundles.
	for(hsize_t i=0; i<bp.second; i++){

		// skip check.
		if( active.find(bp.first[i].ctg_a_idx) == active.end() || active.find(bp.first[i].ctg_b_idx) == active.end() ) continue;

		// get nodes.
		ogdf::node n1 = idx2n[bp.first[i].ctg_a_idx];
		ogdf::node n2 = idx2n[bp.first[i].ctg_b_idx];

		// add node.
		ogdf::edge e = G.newEdge(n1, n2);

	}
}



