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

	// load graph information.
	WRITE_OUT("loading graph information\n");
	NodePair np = node_table_load(node_file);
	BundlePair bp = bundle_table_load(bundle_file);

	// build initial node set.
	NodeSet active;
	for(hsize_t i=0; i<np.second; i++){
		active.insert(np.first[i].node_idx);
	}

	// create initial bundle graph.
	WRITE_OUT("creating initial graph\n");
	BundleGraph BG(np, bp, active);

	// create the initial singular decomp.
	Decomposition decomp;
	decomp.root = -1;
	decomp.node_set = active;
	decomp.decomps = vector<Decomposition>();
	decomp.inters = vector<Intersection>();

	// execute zero decomposition.
	zero_decomp(BG, &decomp);

	// verify zero decomposition.
	verify_connected(np, bp, decomp);

	// digg deeper for 1 decomposition.
	for(hsize_t i=0; i<decomp.decomps.size(); i++){

		// sanity bound.
		if( decomp.decomps[i].node_set.size() < 4 ){
			continue;
		}

		// create the subgraph.
		BundleGraph subg0(np, bp, decomp.decomps[i].node_set);

		// execute one decomposition.
		one_decomp(subg0, &decomp.decomps[i]);

		// verify one decomposition.
		verify_connected(np, bp, decomp.decomps[i]);

		// digg deeper for 2 decomposition.
		for(hsize_t j=0; j<decomp.decomps[i].decomps.size(); j++){

				// sanity bound.
				if( decomp.decomps[i].decomps[j].node_set.size() < 4 ){
					continue;
				}

				// create the subgraph.
				BundleGraph subg1(np, bp, decomp.decomps[i].decomps[j].node_set);

				// execute one decomposition.
				two_decomp(subg1, &decomp.decomps[i].decomps[j]);

				// verify one decomposition.
				verify_connected(np, bp, decomp.decomps[i].decomps[j]);

			}


	}


	return EXIT_SUCCESS;
}

// zero decomposition.
void zero_decomp(BundleGraph BG, Decomposition * decomp){

	// announce.
	WRITE_OUT("executing zero decomposition\n");

	// get graph.
	const Graph G = BG.G;

	// find connected components.
	NodeArray<int> zero_comps(G);
	int cnt = connectedComponents(G, zero_comps);

	// create vector of sets.
	for(int i=0; i<cnt; i++){

		// create decompsition.
		Decomposition tmp;
		tmp.node_set = NodeSet();
		tmp.decomps = vector<Decomposition>();
		tmp.root = -1;

		// add to vector.
		decomp->decomps.push_back(tmp);
	}

	// populate node set.
	node n;
	forall_nodes(n, G){
		decomp->decomps[zero_comps[n]].node_set.insert(BG.n2idx[n]);
	}

}

// one decomposition.
void one_decomp(BundleGraph BG, Decomposition * decomp){

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

	// create vector of sets.
	for(int i=0; i<cnt; i++){
		// create decompsition.
		Decomposition tmp;
		tmp.node_set = NodeSet();
		tmp.decomps = vector<Decomposition>();
		tmp.inters = vector<Intersection>();
		tmp.root = -1;

		// add to vector.
		decomp->decomps.push_back(tmp);
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
				decomp->decomps[cidx].node_set.insert(id1);
				decomp->decomps[cidx].node_set.insert(id2);
			}

			// increment pointer.
			cidx++;
		}
	}

	// verify 1-cuts.
	NodeSet::iterator it;
	vector<int> intersection;
	for(hsize_t i=0; i<decomp->decomps.size(); i++){
		for(hsize_t j=i+1; j<decomp->decomps.size(); j++){

			// clear intersection.
			intersection.clear();

			// check intersection.
			for(it = decomp->decomps[i].node_set.begin(); it!=decomp->decomps[i].node_set.end(); it++){
				if( decomp->decomps[j].node_set.find(*it) != decomp->decomps[j].node_set.end() ){
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

				// create intersection.
				Intersection tmp;
				tmp.compa = i;
				tmp.compb = j;
				tmp.cuts.push_back(intersection.front());

				// insert into vector.
				decomp->inters.push_back(tmp);
			}
		}
	}
}

// two decomposition.
void two_decomp(BundleGraph BG, Decomposition * decomp){

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

	// create vector of sets.
	for(int i=0; i<cnt; i++){
		// create decompsition.
		Decomposition tmp;
		tmp.node_set = NodeSet();
		tmp.decomps = vector<Decomposition>();

		// add to vector.
		decomp->decomps.push_back(tmp);
	}

	// add nodes to sets.
	idx = 0;
	forall_nodes(n, SPQR_T){

		// check for root.
		if( SPQR.rootNode() == n ){
			decomp->root = idx;
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
			decomp->decomps[idx].node_set.insert(id1);
		}

		// increment index.
		idx++;
	}

	// verify 2-cuts.
	NodeSet::iterator it;
	vector<int> intersection;
	for(hsize_t i=0; i<decomp->decomps.size(); i++){
		for(hsize_t j=i+1; j<decomp->decomps.size(); j++){

			// clear intersection.
			intersection.clear();

			// check intersection.
			for(it = decomp->decomps[i].node_set.begin(); it!=decomp->decomps[i].node_set.end(); it++){
				if( decomp->decomps[j].node_set.find(*it) != decomp->decomps[j].node_set.end() ){
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

				// create intersection.
				Intersection tmp;
				tmp.compa = i;
				tmp.compb = j;
				tmp.cuts.push_back(intersection.front());
				tmp.cuts.push_back(intersection.back());

				// insert into vector.
				decomp->inters.push_back(tmp);
			}
		}
	}

}


// verify the decomposition yields connected components.
void verify_connected(NodePair np, BundlePair bp, Decomposition decomp){

	// announce.
	WRITE_OUT("verifying...\n");

	// check its not empty.
	if( decomp.node_set.size() == 0 ){
		WRITE_ERR("error: total decomp is empty.\n");
		exit(1);
	}

	// check each component.
	for(hsize_t i=0; i<decomp.decomps.size(); i++){

		// check its not empty.
		if( decomp.decomps[i].node_set.size() == 0 ){
			WRITE_ERR("error: decomp is empty.\n");
			exit(1);
		}

		// induce subgraph.
		BundleGraph subg(np, bp, decomp.decomps[i].node_set);

		// verify its connected.
		if( isConnected(subg.G) == false ){
			WRITE_ERR("error: is not connected.\n");
			exit(1);
		}
	}
}
