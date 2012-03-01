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

// create global node index.
int GLOBAL_IDX = 0;

int main(int argc, char* argv[]) {

	// arguments.
	const char * node_file = argv[1];
	const char * bundle_file = argv[2];
	unsigned int cutoff = atoi(argv[3]);

	// sanity check.
	if( cutoff < 5 ){
		cutoff = 5;
	}

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

	// verify decomposition.
	verify_connected(np, bp, DG);

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

	}

	// verify decomposition.
	verify_connected(np, bp, DG);

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

	}

	// verify decomposition.
	verify_connected(np, bp, DG);

	// write out decomposition to disk.
	WRITE_OUT("completed succesfully.\n");
	return EXIT_SUCCESS;
}

// zero decomposition.
void zero_decomp(BundleGraph BG, DecompGraph & DG){

	// announce.
	WRITE_OUT("executing zero decomposition\n");

	// get graph.
	const Graph G = BG.G;

	// find connected components.
	NodeArray<int> comps(G);
	int cnt = connectedComponents(G, comps);

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
		DG[vID].global = GLOBAL_IDX;
		GLOBAL_IDX++;
	}

	// populate node set.
	node n;
	forall_nodes(n, G){
		DG[nlookup[comps[n]]].set.insert(BG.n2idx[n]);
	}

}

// one decomposition.
void one_decomp(BundleGraph BG, DecompGraph & DG, VertexID parent){

	// announce.
	WRITE_OUT("executing one decomposition\n");

	// save graph.
	Graph G = BG.G;

	// sanity check.
	if( isConnected(G) == false ){
		WRITE_ERR("not connected: call jwow!\n");
		exit(1);
	}

	// declare vars.
	node n, n1, n2;
	long id1, id2, cidx;

	// find connected components.
	EdgeArray<int> comps(G);
	int cnt = biconnectedComponents(G, comps);

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
		DG[vID].global = GLOBAL_IDX;
		GLOBAL_IDX++;
	}

	// populate node set.
	edge e;
	forall_edges(e, G){

		// get nodes.
		n1 = e->source();
		n2 = e->target();

		// add to graph.
		DG[nlookup[comps[e]]].set.insert(BG.n2idx[n1]);
		DG[nlookup[comps[e]]].set.insert(BG.n2idx[n2]);
	}


	// add nodes to TGraph
	TGraph TG;
	sparse_hash_map<int, VertexID> tlookup;
	for(int i=0; i<cnt; i++){

		// add node.
		TVertexID vID = boost::add_vertex(TG);
		tlookup[i] = vID;
	}

	// add edges to TGraph
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

				// get TGraph nodes.
				TVertexID t1 = tlookup[i];
				TVertexID t2 = tlookup[j];

				// add edge.
				TEdgeID edge;
				bool ok;
				boost::tie(edge, ok) = boost::add_edge(t2, t1, TG);
				if(ok == false){
					WRITE_ERR("problem adding edge.\n");
					exit(1);
				}
				TG[edge].cuts.push_back(intersection.front());
			}
		}
	}

	// DFS this graph to build the decomposition tree.
	MyVisitor vis;
	boost::depth_first_search(TG, boost::visitor(vis));

	// apply the decomposition edges to the actual graph.
	for(vector<Etmp>::iterator it=vis.elist.begin(); it!=vis.elist.end(); it++){

		// grab nodes.
		VertexID p = nlookup[it->source];
		VertexID q = nlookup[it->target];

		// add edges.
		EdgeID edge;
		bool ok;
		boost::tie(edge, ok) = boost::add_edge(p, q, DG);
		if(ok == false){
			WRITE_ERR("problem adding edge.\n");
			exit(1);
		}
		DG[edge].cuts = it->cuts;
	}

	// add root edge.
	TVertexID r = *vertices(TG).first;
	EdgeID edge;
	bool ok;
	boost::tie(edge, ok) = boost::add_edge(parent, nlookup[TG[r].idx], DG);
	if(ok == false){
		WRITE_ERR("problem adding edge.\n");
		exit(1);
	}
}

// two decomposition.
void two_decomp(BundleGraph BG, DecompGraph & DG, VertexID parent){

	// announce.
	WRITE_OUT("executing two decomposition\n");

	// save graph.
	Graph G = BG.G;

	// sanitfy check.
	if( isBiconnected(G) == false ){
		WRITE_ERR("not biconnected: call snooki!\n");
		exit(1);
	}

	// declare vars.
	node n, p, q;
	long idx, id1, id2;

	// get spqr structure.
	StaticSPQRTree SPQR(G);
	const Graph &SPQR_T = SPQR.tree();

	// identify cuts.
	vector<pair<int, int> > cuts;
	forall_nodes(n, SPQR_T){

		// get skeleton and its graph.
		const Skeleton &sk = SPQR.skeleton(n);
		const Graph &g = sk.getGraph();

		// find virtual edges.
		edge e;
		forall_edges(e, g){
			if( sk.isVirtual(e) == true ){
				// get cuts.
				p = sk.original(e->source());
				q = sk.original(e->target());
				id1 = BG.n2idx[p];
				id2 = BG.n2idx[q];

				// add cuts to vector.
				if( id1 < id2 ){
					cuts.push_back(pair<int,int>(id1, id2));
				} else {
					cuts.push_back(pair<int,int>(id2, id1));
				}
			}
		}
	}

	// build list of neighbors for each cut.
	sparse_hash_map<int, vector<int> > neighbors;
	for(hsize_t i=0; i<cuts.size(); i++){

		// simplify.
		id1 = cuts[i].first;
		id2 = cuts[i].second;

		// check if need neighbor search.
		if( neighbors.find(id1) == neighbors.end() ){
			// look for neighbors.
			neighbors[id1] = vector<int>();
			adjEntry adj;
			forall_adj(adj, BG.idx2n[id1]){
				neighbors[id1].push_back(BG.n2idx[adj->theNode()]);
			}
		}

		// check if need neighbor search.
		if( neighbors.find(id2) == neighbors.end() ){
			// look for neighbors.
			neighbors[id2] = vector<int>();
			adjEntry adj;
			forall_adj(adj, BG.idx2n[id2]){
				neighbors[id2].push_back(BG.n2idx[adj->theNode()]);
			}
		}
	}

	// remove all cuts from graph.
	for(sparse_hash_map<int, vector<int> >::iterator it=neighbors.begin(); it!=neighbors.end(); it++){
		G.delNode(BG.idx2n[it->first]);
	}

	// identify connected components.
	NodeArray<int> comps(G);
	int cnt = connectedComponents(G, comps);

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
		DG[vID].global = GLOBAL_IDX;
		GLOBAL_IDX++;
	}

	// populate node set with non-cut.
	sparse_hash_map<int, VertexID> xlookup;
	forall_nodes(n, G){
		DG[nlookup[comps[n]]].set.insert(BG.n2idx[n]);
		xlookup[BG.n2idx[n]] = nlookup[comps[n]];
	}

	// add the cut nodes.
	for(sparse_hash_map<int, vector<int> >::iterator it=neighbors.begin(); it!=neighbors.end(); it++){

		// get cut node.
		id1 = it->first;

		// loop over neighbors.
		for(vector<int>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++){

			// get placed node.
			id2 = *it2;

			// add cut node to placed nodes component.
			DG[xlookup[id2]].set.insert(id1);

		}
	}

	// add nodes to TGraph
	TGraph TG;
	sparse_hash_map<int, VertexID> tlookup;
	for(int i=0; i<cnt; i++){

		// add node.
		TVertexID vID = boost::add_vertex(TG);
		tlookup[i] = vID;
	}

	// add edges to TGraph
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
			if( intersection.size() > 2 ){
				WRITE_ERR("error in two decomposition\n");
				exit(1);
			}

			// sanity check.
			if( intersection.size() == 1 ){
				WRITE_ERR("error in two decomposition2\n");
				exit(1);
			}

			// save intersection.
			if( intersection.size() == 2 ){

				// get TGraph nodes.
				TVertexID t1 = tlookup[i];
				TVertexID t2 = tlookup[j];

				// add edge.
				TEdgeID edge;
				bool ok;
				boost::tie(edge, ok) = boost::add_edge(t2, t1, TG);
				if(ok == false){
					WRITE_ERR("problem adding edge.\n");
					exit(1);
				}
				TG[edge].cuts.push_back(intersection.front());
				TG[edge].cuts.push_back(intersection.back());
			}
		}
	}

	// DFS this graph to build the decomposition tree.
	MyVisitor vis;
	boost::depth_first_search(TG, boost::visitor(vis));

	// apply the decomposition edges to the actual graph.
	for(vector<Etmp>::iterator it=vis.elist.begin(); it!=vis.elist.end(); it++){

		// grab nodes.
		VertexID p = nlookup[it->source];
		VertexID q = nlookup[it->target];

		// add edges.
		EdgeID edge;
		bool ok;
		boost::tie(edge, ok) = boost::add_edge(p, q, DG);
		if(ok == false){
			WRITE_ERR("problem adding edge.\n");
			exit(1);
		}
		DG[edge].cuts = it->cuts;
	}

	// add root edge.
	TVertexID r = *vertices(TG).first;
	EdgeID edge;
	bool ok;
	boost::tie(edge, ok) = boost::add_edge(parent, nlookup[TG[r].idx], DG);
	if(ok == false){
		WRITE_ERR("problem adding edge.\n");
		exit(1);
	}
}

// verify the decomposition yields connected components.
void verify_connected(NodePair np, BundlePair bp, DecompGraph DG){

	// announce.
	WRITE_OUT("verifying...\n");

	int type;
	DecompGraph::vertex_iterator vertexIt, vertexEnd;
	boost::tie(vertexIt, vertexEnd) = vertices(DG);
	for (; vertexIt != vertexEnd; ++vertexIt){

		// dereference vertexIt, get the ID
	    VertexID vertexID = *vertexIt;

	    // get the actual vertex.
	    DGVertex & vertex = DG[vertexID];
	    type = vertex.stage;

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

		// check type.
		if( type > 0 ){
			if( isBiconnected(subg.G) == false ){
				WRITE_ERR("error: is not connected.\n");
				exit(1);
			}
		}
		if( type > 1 ){
			if( isTriconnected(subg.G) == false ){
				WRITE_ERR("error: is not connected.\n");
				exit(1);
			}
		}
	}
}

