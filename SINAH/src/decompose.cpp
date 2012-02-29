/*
 * decompose.cpp
 *
 *  Created on: Feb 27, 2012
 *      Author: jlindsay
 */

// header
#include "decompose.h"

// namespaces
using namespace std;
using namespace ogdf;
using google::sparse_hash_map;
using google::sparse_hash_set;

// declare space for debug text.
char * DEBUGTXT = new char[5012];

int main(int argc, char* argv[]) {

	// arguments.
	const char * node_file = argv[1];
	const char * bundle_file = argv[2];
	unsigned int cutoff = atoi(argv[3]);

	// load graph information.
	WRITE_OUT("loading graph information\n");
	NodePair np = node_table_load(node_file);
	BundlePair bp = bundle_table_load(bundle_file);

	// create the decomposition object.
	DecompGraph DG;

	// build initial node set.
	NodeSet active;
	for(hsize_t i=0; i<np.second; i++){
		active.insert(np.first[i].node_idx);
	}

	// create initial bundle graph.
	WRITE_OUT("creating initial graph\n");
	BundleGraph BG(np, bp, active);

	// execute zero decomposition.
	zero_decomp(BG, DG);

	// verify zero decomposition.
	//verify_connected(np, bp, DG);


	// digg deeper for 1 decomposition.
	DecompGraph::vertex_iterator vertexIt, vertexEnd;
	boost::tie(vertexIt, vertexEnd) = vertices(DG);
	for (; vertexIt != vertexEnd; ++vertexIt) {

		// dereference vertexIt, get the ID
	    VertexID vertexID = *vertexIt;

	    // get the actual vertex.
	    DGVertex & vertex1 = DG[vertexID];

	    // skip non-level 0.
	    if( vertex1.stage != 0 ) continue;

		// sanity bound.
		if( vertex1.set.size() < cutoff ) continue;

		// create the subgraph.
		BundleGraph subg(np, bp, vertex1.set);

		// execute one decomposition.
		one_decomp(subg, DG, vertexID);


		// verify one decomposition.
		//verify_connected(np, bp, DG);
	}

	// digg deeper for 2 decomposition.
	boost::tie(vertexIt, vertexEnd) = vertices(DG);
	for (; vertexIt != vertexEnd; ++vertexIt) {

		// dereference vertexIt, get the ID
	    VertexID vertexID = *vertexIt;

	    // get the actual vertex.
	    DGVertex & vertex1 = DG[vertexID];

	    // skip non-level 1.
	    if( vertex1.stage != 1 ) continue;

		// sanity bound.
		if( vertex1.set.size() < cutoff ) continue;

		// create the subgraph.
		BundleGraph subg(np, bp, vertex1.set);

		// execute one decomposition.
		two_decomp(subg, DG, vertexID);

		// verify one decomposition.
		//verify_connected(np, bp, DG);
	}

	// write out decomposition to disk.

	return EXIT_SUCCESS;
}

// zero decomposition.
void zero_decomp(BundleGraph BG, DecompGraph & DG){

	// announce.
	WRITE_OUT("executing zero decomposition\n");

	// get graph.
	const Graph G = BG.G;

	// find connected components.
	NodeArray<int> zero_comps(G);
	int cnt = connectedComponents(G, zero_comps);

	// add nodes to DG.
	sparse_hash_map<int, VertexID> nlookup;
	for(int i=0; i<cnt; i++){

		// add node.
		VertexID vID = boost::add_vertex(DG);
		nlookup[i] = vID;

		// set properties.
		DG[vID].idx = i;
		DG[vID].stage = 0;
		DG[vID].set.clear();
	}

	// populate node set.
	node n;
	forall_nodes(n, G){
		DG[nlookup[zero_comps[n]]].set.insert(BG.n2idx[n]);
	}

}

// one decomposition.
void one_decomp(BundleGraph BG, DecompGraph & DG, VertexID parent){

	// announce.
	WRITE_OUT("executing one decomposition\n");

	// save graph.
	Graph G = BG.G;

	// declare vars.
	node n, n1, n2;
	long id1, id2, cidx;

	// build BC tree and auxillary graph.
	BCTree BC(G);
	const Graph &BC_T = BC.bcTree();
	const Graph &BC_G = BC.auxiliaryGraph();
	int cnt = BC.numberOfBComps();

	// add nodes to DG.
	sparse_hash_map<int, VertexID> nlookup;
	for(int i=0; i<cnt; i++){

		// add node.
		VertexID vID = boost::add_vertex(DG);
		nlookup[i] = vID;

		// set properties.
		DG[vID].idx = i;
		DG[vID].stage = 1;
		DG[vID].set.clear();

	}

	// create node lookup node array.
	NodeArray<node> node_map(BC_G);

	// assign each Bnode to a component.
	NodeArray<int> bctidx(BC_T);
	cidx = 0;
	forall_nodes(n, BC_T){
		if(BC.typeOfBNode(n) == BC.BComp){

			// map node to component index.
			bctidx[n] = cidx;

			// retrieve edges from b components.
			const SList<edge> x = BC.hEdges(n);

			// add nodes to this component via its edges.
			for(SListConstIterator<edge> it = x.begin(); it!=x.end(); it++){

				// simplify.
				n1 = (*it)->source();
				n2 = (*it)->target();
				id1 = BG.n2idx[BC.original(n1)];
				id2 = BG.n2idx[BC.original(n2)];

				// add nodes to set.
				DG[nlookup[cidx]].set.insert(id1);
				DG[nlookup[cidx]].set.insert(id2);
			}

			// increment pointer.
			cidx++;
		}
	}

	// add root.
	EdgeID edge;
	bool ok;
	boost::tie(edge, ok) = boost::add_edge(parent, nlookup[0], DG);
	if(ok == false){
		WRITE_ERR("problem adding edge.\n");
		exit(1);
	}

	// verify 1-cuts.
	NodeSet::iterator it;
	vector<int> intersection;
	for(int i=0; i<cnt; i++){
		for(int j=i+1; j<cnt; j++){

			// dereference vertexIt, get the ID
			VertexID vertexID1 = nlookup[i];
			VertexID vertexID2 = nlookup[j];

			// get the actual vertex.
			DGVertex & vertex1 = DG[vertexID1];
			DGVertex & vertex2 = DG[vertexID2];

			// clear intersection.
			intersection.clear();

			// check intersection.
			for(it = vertex1.set.begin(); it!=vertex1.set.end(); it++){
				if( vertex2.set.find(*it) != vertex2.set.end() ){
					intersection.push_back(*it);
				}
			}

			// sanity check.
			if( intersection.size() > 1 ){
				WRITE_ERR("error in one decomposition\n");
				exit(1);
			}

			// save intersection.
			if( intersection.size() == 1 ){

				// insert into vector.
				EdgeID edge;
				bool ok;
				boost::tie(edge, ok) = boost::add_edge(vertexID1, vertexID2, DG);
				if(ok == false){
					WRITE_ERR("problem adding edge.\n");
					exit(1);
				}
				DG[edge].cuts.push_back(intersection.front());
			}
		}
	}
}

// two decomposition.
void two_decomp(BundleGraph BG, DecompGraph & DG, VertexID parent){

	// announce.
	WRITE_OUT("executing two decomposition\n");

	// save graph.
	Graph G = BG.G;

	// declare vars.
	node n, p, q;
	long idx, id1;

	// get spqr structure.
	StaticSPQRTree SPQR(G);
	const Graph &SPQR_T = SPQR.tree();
	int cnt = SPQR_T.numberOfNodes();

	// add nodes to DG.
	sparse_hash_map<int, VertexID> nlookup;
	for(int i=0; i<cnt; i++){

		// add node.
		VertexID vID = boost::add_vertex(DG);
		nlookup[i] = vID;

		// set properties.
		DG[vID].idx = i;
		DG[vID].stage = 2;
		DG[vID].set.clear();

	}

	// add nodes to sets.
	idx = 0;
	forall_nodes(n, SPQR_T){

		// check for root.
		if( SPQR.rootNode() == n ){
			EdgeID edge;
			bool ok;
			boost::tie(edge, ok) = boost::add_edge(parent, nlookup[idx], DG);
			if(ok == false){
				WRITE_ERR("problem adding edge.\n");
				exit(1);
			}
		}

		// get skeleton and its graph.
		const Skeleton &sk = SPQR.skeleton(n);
		const Graph &g = sk.getGraph();

		// verify skeleton is connected.
		if( isConnected(g) == false ){
			printf("error skeleton is not connected.\n");
			exit(1);
		}

		// add nodes to set.
		forall_nodes(p, g){
			q = sk.original(p);
			id1 = BG.n2idx[q];
			DG[nlookup[idx]].set.insert(id1);
		}

		// increment index.
		idx++;
	}

	// verify 2-cuts.
	NodeSet::iterator it;
	vector<int> intersection;
	for(int i=0; i<cnt; i++){
		for(int j=i+1; j<cnt; j++){

			// dereference vertexIt, get the ID
			VertexID vertexID1 = nlookup[i];
			VertexID vertexID2 = nlookup[j];

			// get the actual vertex.
			DGVertex & vertex1 = DG[vertexID1];
			DGVertex & vertex2 = DG[vertexID2];

			// clear intersection.
			intersection.clear();

			// check intersection.
			for(it = vertex1.set.begin(); it!=vertex1.set.end(); it++){
				if( vertex2.set.find(*it) != vertex2.set.end() ){
					intersection.push_back(*it);
				}
			}

			// sanity check.
			if( intersection.size() > 2  ){
				WRITE_ERR("error in two decomposition\n");
				for(vector<int>::iterator it2 = intersection.begin(); it2!= intersection.end(); it2++){
					cout << *it2 << endl;
				}
				exit(1);
			}

			// save intersection.
			if( intersection.size() == 2 ){

				// insert into vector.
				EdgeID edge;
				bool ok;
				boost::tie(edge, ok) = boost::add_edge(vertexID1, vertexID2, DG);
				if(ok == false){
					WRITE_ERR("problem adding edge.\n");
					exit(1);
				}
				DG[edge].cuts.push_back(intersection.front());
				DG[edge].cuts.push_back(intersection.back());

			}
		}
	}

}

// verify the decomposition yields connected components.
void verify_connected(NodePair np, BundlePair bp, DecompGraph DG){

	// announce.
	WRITE_OUT("verifying...\n");

	DecompGraph::vertex_iterator vertexIt, vertexEnd;
	boost::tie(vertexIt, vertexEnd) = vertices(DG);
	for (; vertexIt != vertexEnd; ++vertexIt){

		// dereference vertexIt, get the ID
	    VertexID vertexID = *vertexIt;

	    // get the actual vertex.
	    DGVertex & vertex = DG[vertexID];

	    // check if empty.
	    if( vertex.set.size() == 0 ){
			WRITE_ERR("error: decomp is empty.\n");
			exit(1);
	    }

		// induce subgraph.
		BundleGraph subg(np, bp, vertex.set);

		// verify its connected.
		if( isConnected(subg.G) == false ){
			WRITE_ERR("error: is not connected.\n");
			exit(1);
		}
	}
}

