/*
 * decompose.h
 *
 *  Created on: Feb 27, 2012
 *      Author: jlindsay
 */

#ifndef DECOMPOSE_H_
#define DECOMPOSE_H_

// google
#include <google/sparse_hash_map>
#include <google/sparse_hash_set>

// hdf5
#include "types.h"

// bundle graph.
#include "BundleGraph.h"

// logging.
#include "logging.h"

// used to store intersections between decompositions.
typedef struct Intersection {

	// the type of intersection: 1, 2, 3+
	int type;

	// the two components involved in intersection.
	int compa;
	int compb;

	// list of cut vertex.
	std::vector<int> cuts;

} Intersection;

// used to store the decomposition.
typedef struct Decomposition {

	// the active node set.
	NodeSet node_set;

	// a vector of subsets.
	std::vector<Decomposition> decomps;

	// a vector of intersections.
	std::vector<Intersection> inters;

	// root of intersection if it can be rooted.
	int root;
} Decomposition;

// decomposition functions.
void zero_decomp(BundleGraph BG, Decomposition * decomp);
void one_decomp(BundleGraph BG, Decomposition * decomp);
//void two_decomp(BundleGraph BG, Decomposition * decomp);
void two_decomp(BundleGraph BG, NodePair np, BundlePair bp, Decomposition * decomp);

// ancillary functions.
void verify_connected(NodePair np, BundlePair bp, Decomposition decomp);

#endif /* DECOMPOSE_H_ */
