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

// boost BGL
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>

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

// used to store a result for return.
typedef struct Result {
	std::vector<Decomposition> decomps;
	std::vector<Intersection> inters;
	int root;
} Result;

struct DGVertex {
    int idx;
	int gidx;
	int nidx;
};

struct DGEdge {
    std::vector<int> cuts;
};

typedef boost::adjacency_list<  // adjacency_list is a template depending on :
    boost::listS,               //  The container used for egdes : here, std::list.
    boost::vecS,                //  The container used for vertices: here, std::vector.
    boost::undirectedS,           //  directed or undirected edges ?.
    DGVertex,                     //  The type that describes a Vertex.
    DGEdge                        //  The type that describes an Edge
> DecompGraph;

typedef DecompGraph::vertex_descriptor VertexID;
typedef DecompGraph::edge_descriptor   EdgeID;

// decomposition functions.
void zero_decomp(BundleGraph BG, DecompGraph & DG, std::vector<DecompGraph> & VDG, std::vector<NodeSet> & VNS);
void one_decomp(BundleGraph BG, DecompGraph & DG, std::vector<DecompGraph> & VDG, std::vector<NodeSet> & VNS);
void two_decomp(BundleGraph BG, DecompGraph & DG, std::vector<DecompGraph> & VDG, std::vector<NodeSet> & VNS);

// ancillary functions.
void verify_connected(NodePair np, BundlePair bp, DecompGraph decomp, std::vector<NodeSet> VNS);

#endif /* DECOMPOSE_H_ */
