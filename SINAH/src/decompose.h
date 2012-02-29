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

struct DGVertex {
    int idx;
    int stage;
	google::sparse_hash_set<int> set;
};

struct DGEdge {
    std::vector<int> cuts;
};

typedef boost::adjacency_list<  // adjacency_list is a template depending on :
    boost::listS,               //  The container used for egdes : here, std::list.
    boost::vecS,                //  The container used for vertices: here, std::vector.
    boost::directedS,           //  directed or undirected edges ?.
    DGVertex,                     //  The type that describes a Vertex.
    DGEdge                        //  The type that describes an Edge
> DecompGraph;

typedef DecompGraph::vertex_descriptor VertexID;
typedef DecompGraph::edge_descriptor   EdgeID;

// decomposition functions.
void zero_decomp(BundleGraph BG, DecompGraph & DG);
void one_decomp(BundleGraph BG, DecompGraph & DG, VertexID parent);
void two_decomp(BundleGraph BG, DecompGraph & DG, VertexID parent);

// ancillary functions.
void verify_connected(NodePair np, BundlePair bp, DecompGraph decomp);

#endif /* DECOMPOSE_H_ */
