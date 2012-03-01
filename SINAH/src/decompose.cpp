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
	const char * result_dir = argv[3];
	unsigned int cutoff = atoi(argv[4]);

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

	// build vector of vertexID's.
	vector<VertexID> buffer;
	DecompGraph::vertex_iterator vertexIt, vertexEnd;
	boost::tie(vertexIt, vertexEnd) = vertices(DG);
	for (; vertexIt != vertexEnd; ++vertexIt) {
		buffer.push_back(*vertexIt);
	}

	// digg deeper for 1 decomposition.
	for(unsigned int i=0; i<buffer.size(); i++){

		// dereference vertexIt, get the ID
	    VertexID vertexID = buffer[i];

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

	// build vector of vertexID's.
	buffer.clear();
	boost::tie(vertexIt, vertexEnd) = vertices(DG);
	for (; vertexIt != vertexEnd; ++vertexIt) {
		buffer.push_back(*vertexIt);
	}

	// digg deeper for 2 decomposition.
	for(unsigned int i=0; i<buffer.size(); i++){

		// dereference vertexIt, get the ID
	    VertexID vertexID = buffer[i];

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

    // make the annotation files.
    char zero_file[500];
    char one_file[500];
    char two_file[500];

    sprintf(zero_file, "%s/decomp_zero.txt", result_dir);
    sprintf(one_file, "%s/decomp_one.txt", result_dir);
    sprintf(two_file, "%s/decomp_two.txt", result_dir);

    // open the annotation files.
    FILE * zero_out = fopen(zero_file, "w");
    FILE * one_out = fopen(one_file, "w");
    FILE * two_out = fopen(two_file, "w");

	// build vector of vertexID's.
	buffer.clear();
	boost::tie(vertexIt, vertexEnd) = vertices(DG);
	for (; vertexIt != vertexEnd; ++vertexIt) {
		buffer.push_back(*vertexIt);
	}

	// build vector of edgeID's
	vector<EdgeID> ebuffer;
	DecompGraph::edge_iterator edgeIt, edgeEnd;
	BOOST_FOREACH(EdgeID eid, edges(DG)) {
		ebuffer.push_back(eid);
	}

	// build parent lookup.
	sparse_hash_map<int, VertexID> glookup;
	for(unsigned int i=0; i<buffer.size(); i++){
		glookup[DG[buffer[i]].global] = buffer[i];
	}

	// write zero.
	for(unsigned int i=0; i<buffer.size(); i++){

		// prep.
	    VertexID vertexID = buffer[i];
	    DGVertex & vertex1 = DG[vertexID];

	    // skip non-level 0.
	    if( vertex1.stage != 0 ) continue;

		// write each zero node.
		fprintf(zero_out, "stage=0\ttype=N\tidx=%d,,\t", vertex1.idx);

		// fill the set.
		fprintf(zero_out, "set=");
		for(NodeSet::iterator it=vertex1.set.begin(); it!=vertex1.set.end(); it++){
			fprintf(zero_out, "%d,", *it);
		}
		fprintf(zero_out, "\n");
	}

	// write one.
	int ida1, ida2, ida3, idb1, idb2, idb3;
	for(unsigned int i=0; i<buffer.size(); i++){

		// prep.
	    VertexID vertexID = buffer[i];
	    DGVertex & vertex1 = DG[vertexID];

	    // skip non-level 1.
	    if( vertex1.stage != 1 ) continue;

	    // prep 2.
	    ida1 = DG[glookup[vertex1.pidx]].idx;
	    ida2 = vertex1.idx;

		// write each one node.
		fprintf(one_out, "stage=1\ttype=N\tidx=%d,%d\t", ida1, ida2);

		// fill the set.
		fprintf(one_out, "set=");
		for(NodeSet::iterator it=vertex1.set.begin(); it!=vertex1.set.end(); it++){
			fprintf(one_out, "%d,", *it);
		}
		fprintf(one_out, "\n");
    }

    // write one edges.
	int pid1, pid2, c1;
	for(unsigned int i=0; i<ebuffer.size(); i++){

		// get source and target info.
		pid1 = DG[source(ebuffer[i], DG)].pidx;
		pid2 = DG[target(ebuffer[i], DG)].pidx;

		// only add edges in same component.
		if( pid1 != pid2 ) continue;

		// only add edges in 1 stage.
		if( DG[target(ebuffer[i], DG)].stage != 1 ) continue;

		// get parents local idx.
		ida1 = DG[glookup[pid1]].idx;
		ida2 = DG[source(ebuffer[i], DG)].idx;

		idb1 = DG[glookup[pid2]].idx;
		idb2 = DG[target(ebuffer[i], DG)].idx;

		// get singular cut.
		c1 = DG[ebuffer[i]].cuts.front();

		// write each one edge.
		fprintf(one_out, "stage=1\ttype=E\tidx1=%d,%d\tidx2=%d,%d\tcut=%d\n", ida1, ida2, idb1, idb2, c1);
    }

	// write two.
	int id3;
	for(unsigned int i=0; i<buffer.size(); i++){

		// prep.
	    VertexID vertexID = buffer[i];
	    DGVertex & vertex1 = DG[vertexID];

	    // skip non-level 2.
	    if( vertex1.stage != 2 ) continue;

	    // prep 2.
	    ida1 = DG[DG[glookup[vertex1.pidx]].pidx].idx;
	    ida2 = DG[glookup[vertex1.pidx]].idx;
	    ida3 = vertex1.idx;

		// write each one node.
		fprintf(two_out, "stage=2\ttype=N\tidx=%d,%d,%d\t", ida1, ida2, ida3);

		// fill the set.
		fprintf(two_out, "set=");
		for(NodeSet::iterator it=vertex1.set.begin(); it!=vertex1.set.end(); it++){
			fprintf(two_out, "%d,", *it);
		}
		fprintf(two_out, "\n");
    }

    // write two edges.
	int pid3, pid4, c2;
	for(unsigned int i=0; i<ebuffer.size(); i++){

		// get source and target info.
		pid1 = DG[source(ebuffer[i], DG)].pidx;
		pid2 = DG[target(ebuffer[i], DG)].pidx;
		pid3 = DG[glookup[pid1]].pidx;
		pid4 = DG[glookup[pid2]].pidx;

		// only add edges in same component.
		if( pid1 != pid2 || pid3 != pid4 ) continue;

		// only add edges in 2 stage.
		if( DG[target(ebuffer[i], DG)].stage != 2 ) continue;

		// get parents local idx.
		ida1 = DG[glookup[DG[glookup[pid1]].pidx]].pidx;
		ida2 = DG[glookup[pid1]].idx;
		ida3 = DG[source(ebuffer[i], DG)].idx;

		idb1 = DG[glookup[DG[glookup[pid2]].pidx]].pidx;
		idb2 = DG[glookup[pid2]].idx;
		idb3 = DG[target(ebuffer[i], DG)].idx;

		// get singular cut.
		c1 = DG[ebuffer[i]].cuts.front();
		c2 = DG[ebuffer[i]].cuts.back();

		// write each one edge.
		fprintf(two_out, "stage=2\ttype=E\tidx1=%d,%d,%d\tidx2=%d,%d,%d\tcut=%d,%d\n", ida1, ida2, ida3, idb1, idb2, idb3, c1, c2);
    }

	// write roots.
	for(unsigned int i=0; i<ebuffer.size(); i++){

		// get source and target info.
		pid1 = DG[source(ebuffer[i], DG)].pidx;
		pid2 = DG[target(ebuffer[i], DG)].pidx;

		// look for target in 2.
		if( DG[target(ebuffer[i], DG)].stage != 2 ) continue;

		// look for different source.
		if( DG[source(ebuffer[i], DG)].stage == 1 ) {

			// target is cut.
			idb1 = DG[glookup[DG[glookup[pid2]].pidx]].pidx;
			idb2 = DG[glookup[pid2]].idx;
			idb3 = DG[target(ebuffer[i], DG)].idx;

			// target is cut.
			 fprintf(two_out, "stage=2\ttype=R\tidx=%d,%d,%d\n", idb1, idb2, idb3);
		}
	}

	// close files.
    fclose(zero_out);
    fclose(one_out);
    fclose(two_out);

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
		DG[vID].pidx = -1;
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
		DG[vID].pidx = DG[parent].global;
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

		// annotate it.
		TG[vID].idx = i;
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
	vector<Etmp> elist;
	MyVisitor vis(&elist);
	boost::depth_first_search(TG, boost::visitor(vis));

	// apply the decomposition edges to the actual graph.
	for(vector<Etmp>::iterator it=elist.begin(); it!=elist.end(); it++){

		// grab nodes.
		VertexID p = nlookup[it->source];
		VertexID q = nlookup[it->target];

		// sanity check.
		if(it->source == it->target){
			WRITE_ERR("BAD DFS1\n");
			exit(1);
		}

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

	// copy graph for messing with.
	GraphCopySimple GC(G);

	// remove all cuts from graph.
	for(sparse_hash_map<int, vector<int> >::iterator it=neighbors.begin(); it!=neighbors.end(); it++){
		GC.delNode(GC.copy(BG.idx2n[it->first]));
	}

	// identify connected components.
	NodeArray<int> comps(GC);
	int cnt = connectedComponents(GC, comps);

	// add nodes to DG.
	sparse_hash_map<int, VertexID> nlookup;
	for(int i=0; i<cnt; i++){

		// add node.
		VertexID vID = boost::add_vertex(DG);
		nlookup[i] = vID;

		// set properties.
		DG[vID].idx = i;
		DG[vID].pidx = DG[parent].global;
		DG[vID].stage = 2;
		DG[vID].set.clear();
		DG[vID].global = GLOBAL_IDX;
		GLOBAL_IDX++;
	}

	// populate node set with non-cut.
	sparse_hash_map<int, VertexID> xlookup;
	forall_nodes(n, GC){
		DG[nlookup[comps[n]]].set.insert(BG.n2idx[GC.original(n)]);
		xlookup[BG.n2idx[GC.original(n)]] = nlookup[comps[n]];
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

		// annotate it.
		TG[vID].idx = i;
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
	vector<Etmp> elist;
	MyVisitor vis(&elist);
	boost::depth_first_search(TG, boost::visitor(vis));

	// apply the decomposition edges to the actual graph.
	for(vector<Etmp>::iterator it=elist.begin(); it!=elist.end(); it++){

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
			WRITE_ERR("error: is not connected0.\n");
			exit(1);
		}

	}
}

