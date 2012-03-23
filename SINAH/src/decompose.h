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
#include <boost/graph/depth_first_search.hpp>
#include <boost/foreach.hpp>
#include <boost/parameter/name.hpp>

// bundle graph.
#include "BundleGraph.h"

// logging.
#include "logging.h"

typedef struct DGVertex {
    int idx;
    int stage;
    int global;
	google::sparse_hash_set<int> set;
} DGVertex;

typedef struct TVertex {
    int idx;
} TVertex;

typedef struct DGEdge {
    std::vector<int> cuts;
} DGEdge;

typedef boost::adjacency_list<  // adjacency_list is a template depending on :
    boost::listS,               //  The container used for egdes : here, std::list.
    boost::vecS,                //  The container used for vertices: here, std::vector.
    boost::directedS,           //  directed or undirected edges ?.
    DGVertex,                     //  The type that describes a Vertex.
    DGEdge                        //  The type that describes an Edge
> DecompGraph;

typedef boost::adjacency_list<  // adjacency_list is a template depending on :
    boost::listS,               //  The container used for egdes : here, std::list.
    boost::vecS,                //  The container used for vertices: here, std::vector.
    boost::undirectedS,           //  directed or undirected edges ?.
    TVertex,                     //  The type that describes a Vertex.
    DGEdge                        //  The type that describes an Edge
> TGraph;

typedef DecompGraph::vertex_descriptor VertexID;
typedef DecompGraph::edge_descriptor   EdgeID;

typedef TGraph::vertex_descriptor TVertexID;
typedef TGraph::edge_descriptor   TEdgeID;

typedef google::sparse_hash_map<int, VertexID> ParentMap;

// decomposition functions.
void zero_decomp(BundleGraph BG, DecompGraph & DG);
void one_decomp(BundleGraph BG, DecompGraph & DG, ParentMap & PMAP, VertexID parent);
void two_decomp(BundleGraph BG, DecompGraph & DG, ParentMap & PMAP, VertexID parent);

// ancillary functions.
void verify_connected(NodePair np, BundlePair bp, DecompGraph decomp);

// visitor for DFS.
typedef struct Etmp {
	int source;
	int target;
    std::vector<int> cuts;
} Etmp;
class MyVisitor : public boost::default_dfs_visitor {
public:

	// constructor.
	MyVisitor(std::vector<Etmp> * e){
		elist = e;
	}

	// track directed edges in tree.
	std::vector<Etmp> * elist;
	void tree_edge(TEdgeID e, const TGraph& g) {

		// get parent and child.
		TVertexID p = source(e, g);
		TVertexID q = target(e, g);

		// build entry.
		Etmp tmp;
		tmp.source = g[p].idx;
		tmp.target = g[q].idx;
		tmp.cuts = g[e].cuts;

		// add it.
		elist->push_back(tmp);

		return;
	}
};

#endif /* DECOMPOSE_H_ */
