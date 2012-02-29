'''
loads all relavent information and solves ILP
'''
# program imports.
from data_structs.nodes import load_nodes
from data_structs.edges import load_edges
from data_structs.bundles import load_bundles
from data_structs.agps import load_agps
from data_structs.types import agp_dt

from graphs.agp import ScaffAgp

# system imports.
import math
import numpy as np
import subprocess
import networkx as nx
import logging
import sys
import os

# logging.
#logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )

# parameters.
test_file = os.path.abspath(sys.argv[1])
ref_file = os.path.abspath(sys.argv[2])
name = sys.argv[3]
method = sys.argv[4]

########### classes ####################

########### functions ##################
	
def make_agp(fpath):
	''' returns agp graph '''
	
	# load agp array.
	agp_edges = load_agps(fpath)
	
	# build AGP graph.
	AG = ScaffAgp(agp_edges)	
	
	return AG, agp_edges

def build_adj(G):
	''' build adjacent set '''
	
	# create empty set.
	adj = set()
	
	# loop over edges and add.
	for e in G.edges():
		
		# add to set.
		adj.add(e)
		
	# return set.
	return adj

def compute_whole(nodes, true_adj, TG):
	''' computes statistics when TG is already oriented. This
	method explicitly represents all 4 sets and should only be
	used to test the logic of these operations or on some
	small test cases to understand error and stuff.'''
	
	# build false adj.
	false_adj = set()
	for n1 in nodes:
		for n2 in nodes:
			
			# skip same and true_adj.
			if n1 == n2 or (n1, n2) in true_adj:
				continue
				
			# add false adj.
			false_adj.add((n1, n2))
	
	# build predicted adj.
	pred_adj = build_adj(TG)
		
	# build predicted nonadj.
	pred_nonadj = set()
	for n1 in nodes:
		for n2 in nodes:
			
			# skip same and true_adj.
			e = (n1, n2)
			if n1 == n2 or e in pred_adj:
				continue
				
			# add false adj.
			pred_nonadj.add(e)

	# compute the A statistic set.
	A = true_adj.intersection(pred_adj)
		
	# compute the B statistic set.
	B = false_adj.intersection(pred_adj)
	
	# compute the C statistic set.
	C = true_adj.intersection(pred_nonadj)

	# compute the D statistic set.
	D = false_adj.intersection(pred_nonadj)
		
	# return categories.
	return len(A), len(B), len(C), len(D)
	
def compute_formula(nodes, true_adj, RG, TG):
	''' computes statistics when TG is already oriented. This
	method computes A and B exactly, but uses formulas
	to find C and D.'''
	
	# set vars.
	n = RG.number_of_nodes()
	m = RG.number_of_edges()
	N = n * (n-1)
	
	# build predicted adj.
	pred_adj = build_adj(TG)

	# compute the A statistic set.
	A = true_adj.intersection(pred_adj)
		
	# compute the B statistic set.
	B = pred_adj.difference(true_adj)
	
	# compute the C statistic set.
	C = m - len(A)

	# compute the D statistic set.
	D = N - m - len(B)
		
	# return categories.
	return len(A), len(B), C, D, A

def compute_orientation(ref_agp, test_agp, TG, true_adj):
	''' calculates orientation based on maximized A '''
	
	# track global good and bad.
	good = 0
	bad = 0
	
	# map contig to reference idx.
	ref_idx = {}
	for i in range(ref_agp.size):
		ref_idx[ref_agp[i]['comp_name']] = i

	# map contig to reference idx.
	test_idx = {}
	for i in range(test_agp.size):
		test_idx[test_agp[i]['comp_name']] = i

	# loop over each component.
	for comp in nx.weakly_connected_components(TG):
		
		# build predicted adjacents (both ways).
		pred_adj_1 = set()
		pred_adj_2 = set()
		for e in TG.edges(comp):
			pred_adj_1.add( (e[0], e[1]) )
			pred_adj_2.add( (e[1], e[0]) )
			
		# compute the A statistic set.
		A_1 = true_adj.intersection(pred_adj_1)
		A_2 = true_adj.intersection(pred_adj_2)

		# test which has best A statistic.
		if len(A_1) < len(A_2):
			
			# flip entries in test AGP array.
			for n in comp:
				test_agp[test_idx[n]]['comp_orien'] = 1 - test_agp[test_idx[n]]['comp_orien']
				
		# compute sizes.
		ref_size = [0, 0]
		test_size = [0, 0]
				
		for n in comp:
			ref_size[ ref_agp[ ref_idx[n] ]['comp_orien'] ] += 1
			test_size[ test_agp[ test_idx[n] ]['comp_orien'] ] += 1
			
		# find the difference.
		tmp1 = sorted(ref_size)
		tmp2 = sorted(test_size)
		
		# note it.
		bad += abs(tmp1[0] - tmp2[0])
		good += len(comp) - abs(tmp1[0] - tmp2[0])
		
					
	# return counts.
	return good, bad
			
def orient_TG(nodes, true_adj, TG):
	''' creates a new TG with edges orientated by max category A'''

	# create new directed graph.
	NG = nx.DiGraph()
	
	# add all nodes.
	for n in nodes:
		NG.add_node(n)

	# loop over components.
	elist = []
	for comp in nx.weakly_connected_components(TG):

		# build predicted adjacents (both ways).
		pred_adj_1 = set()
		pred_adj_2 = set()
		for e in TG.edges(comp):
			pred_adj_1.add( (e[0], e[1]) )
			pred_adj_2.add( (e[1], e[0]) )
			
		# compute the A statistic set.
		A_1 = true_adj.intersection(pred_adj_1)
		A_2 = true_adj.intersection(pred_adj_2)

		# test which has best A statistic.
		if len(A_1) > len(A_2):
			A = A_1
			pred_adj = pred_adj_1
		else:
			A = A_2
			pred_adj = pred_adj_2
			
		# add edges to new graph.
		for e in pred_adj:
			elist.append(e)
			
	# add edges in one fell swoop.
	NG.add_edges_from(elist)
			
	# return new graph.
	return NG
			

	
def sanity_check(RG, TG, A, B, C, D):
	''' make sure these parameter make sense '''
	
	# set vars.
	n = RG.number_of_nodes()
	m = RG.number_of_edges()
	N = n * (n-1)
	
	# run check 1: m = A + C.
	check1 = m == A + C
	
	# run check 2: N - m = B + D.
	check2 = N - m == B + D
	
	# run report.
	fail = False
	
	if check1 == False:
		logging.error("failed check 1: m = A + C, %i = %i + %i: %i" % (m, A, C, m == A + C))
		fail = True

	if check2 == False:
		logging.error("failed check 2: N-m = B + D, %i - %i = %i + %i: %i" % (N, m, B, D, N-m==B+D))
		fail = True
		
	if fail == True:
		sys.exit(1)
		
def ensure_consistent(TG, nodes, ref_edges, test_edges):
	''' makes sure they contain same nodes '''
	
	# build index of nodes.
	idx = {}
	for i in range(ref_edges.size):
		if ref_edges[i]['comp_type'] == "W":
			idx[ref_edges[i]['comp_name']] = i
	
	# add missing nodes.
	for n in nodes:
		
		# node isn't present.
		if TG.has_node(n) == False:
			
			# add it to graph.
			TG.add_node(n)
			
			# resize the numpy array.
			test_edges = np.resize(test_edges, test_edges.size + 1)
			
			# check if there.
			if n not in idx:
				logging.error("wow: n not in idx")
				logging.error("n: %s\n" % (n))
				continue
				#sys.exit(1)
			if idx[n] not in ref_edges:
				logging.error("wow: idx[n] not in ref_edges")
				logging.error("n: %s\tidx[n]:%s\n" % (n, idx[n]))
				continue
				#sys.exit(1)
			
			# copy info into it.
			test_edges[-1] = ref_edges[idx[n]]
			
			# get next scaffold id.
			scaf_name = "scaffold_%i" % (int(test_edges[-2]['scaf_name'].split("_")[1]) + 1)
			
			# update other vars.
			test_edges[-1]['scaf_name'] = scaf_name
			test_edges[-1]['scaf_idx'] = 1
			test_edges[-1]['scaf_start'] = 1
			test_edges[-1]['scaf_stop'] = test_edges[-1]['comp_stop']
			
				
	# make set out of nodes.
	nset = set(nodes)
			
	# check for extra nodes.
	fail = False
	for n in TG.nodes():
		
		# error check.
		if n not in nset:
			fail = True
			logging.error("extra node: %s" % n)
			
	# die.
	if fail == True:
		sys.exit(1)
		
	# return TG and test edges.
	return TG, test_edges
	
def gap_deviation(A_SET, ref_edges, test_edges):
	''' calculates the deviation from the actual gap size for the A
	set of edges '''
	
	# build dictionary for edge distance.
	act_dist = {}
	
	# identify idx of edges in ref.
	for i in range(ref_edges.size - 2):
		
		# check for contig.
		if ref_edges[i]['comp_type'] != "W":
			continue
			
		# check the two adjacent are in the same scaffold.
		if ref_edges[i]['scaf_name'] != ref_edges[i+2]['scaf_name']:
			continue
		
		# check gap between them.
		if ref_edges[i+1]['comp_type'] != "N":
			
			
			logging.error("error in AGP file:")			
			print '\t'.join([str(x) for x in ref_edges[i]])
			print '\t'.join([str(x) for x in ref_edges[i+1]])
			print '\t'.join([str(x) for x in ref_edges[i+2]])
			
			logging.error("*")
			logging.error("*")
			logging.error("*")
			logging.error("*")
			sys.exit(1)
			
		# get dist.
		edist = ref_edges[i+1]['comp_stop'] - ref_edges[i+1]['comp_start']
		if edist == 0:
			edist = 1
			
		# check if this edge is in active set.
		e1 = ( ref_edges[i]['comp_name'], ref_edges[i+2]['comp_name'] )
		e2 = ( ref_edges[i+2]['comp_name'], ref_edges[i]['comp_name'] )

		# add note to appropriate edge.
		if e1 in A_SET:
			act_dist[e1] = edist
		elif e2 in A_SET:
			act_dist[e2] = edist

	# idnetify edges in test set.
	dev = 0.0
	cnt = 0.0
	for i in range(test_edges.size - 2):
		
		# check for contig.
		if test_edges[i]['comp_type'] != "W":
			continue
			
		# check the two adjacent are in the same scaffold.
		if test_edges[i]['scaf_name'] != test_edges[i+2]['scaf_name']:
			continue
		
		# check gap between them.
		if test_edges[i+1]['comp_type'] != "N":
			logging.error("error in AGP file:\n%s\n" % '\t'.join([str(x) for x in test_edges[i+1]]))
			sys.exit(1)
					
		# get dist.
		edist = test_edges[i+1]['comp_stop'] - test_edges[i+1]['comp_start']
			
		# check if this edge is in active set.
		e1 = ( test_edges[i]['comp_name'], test_edges[i+2]['comp_name'] )
		e2 = ( test_edges[i+2]['comp_name'], test_edges[i]['comp_name'] )
		
		if e1 in act_dist:
			adist = act_dist[e1]
		elif e2 in act_dist:
			adist = act_dist[e2]
		else:
			continue
		
		dev += float(abs(adist - edist)) / float(adist)
		cnt += 1.0
		
	# compute average.
	return dev / cnt
		
def get_runtime(agp_file):
	
	# read in agp.
	fin = open(agp_file, "rb")
	lines = fin.readlines()
	fin.close()	
	
	if lines[-1].count("RUNTIME") == 0:
		return -1.0
	
	# parse last line.
	tmp = lines[-1].strip().split()
	
	return float(tmp[1])		
		
def calculate_n50(edges):
	''' calculates scaffold N50 given agp list'''
	
	# count the number of scaffolds.
	ctgset = {}
	ctgsize = 0
	for i in range(edges.size):
		if edges[i]['scaf_name'] not in ctgset:
			ctgset[ edges[i]['scaf_name'] ] = ctgsize
			ctgsize += 1
		
	# create an array of scaffold size integers.
	ctgarr = np.zeros(ctgsize, dtype=np.int)
	for i in range(ctgsize):
		ctgarr[i] = -1
	
	# annotate the array.
	for i in range(edges.size):
		ctgarr[ ctgset[ edges[i]['scaf_name'] ] ] += edges[i]['scaf_stop'] - edges[i]['scaf_start']
	
	# sort the array.
	ctgarr.sort()
	
	# calculate total sequence length.
	totsize = np.sum(ctgarr)
	
	# calculate middle.
	middle = int( float(totsize) / 2.0 )
	
	# find the n50 from the back.
	for i in range(ctgarr.size - 1, -1, -1):
		
		# check if we below 50.
		if totsize <= middle:
			n50 = ctgarr[i]
			break
		
		# make subtraction.
		totsize -= ctgarr[i]
		
	# return N50
	return n50
			
def trueify_scaffold(edges, TG, a_set):
	''' breaks scaffolds at bad edges '''
	
	# count bad edges.
	scafs = list()
	active = set()
	bad = 0
	pscaf = edges[0]['scaf_name']
	for i in range(edges.size):
		
		# check for scaf break.
		if pscaf != edges[i]['scaf_name']:
			scafs.append(active)
			active = set()
			pscaf = edges[i]['scaf_name']
		
		# just add contigs.
		if edges[i]['comp_type'] == "W":
			active.add(i)
			continue
		
		# find edges.
		idx1 = edges[i-1]['comp_name']
		idx2 = edges[i+1]['comp_name']
		e1 = (idx1, idx2)
		e2 = (idx2, idx1)
		
		# add good edges to active.
		if e1 in a_set or e2 in a_set: 
			active.add(i)
			continue
			
		# append active to set.
		scafs.append(active)
		active = set()
		bad += 1
		
		
	# allocate new edge array.
	ne = np.zeros(edges.size - bad, dtype=agp_dt)
	
	# copy data to this.
	idx = 0
	scaf_idx = 0
	for scaf in scafs:
		
		# add these edges.
		for i in scaf:
			
			# add edge.
			ne[idx] = edges[i]
			
			# update name.
			ne[idx]['scaf_name'] = 'scaf_%i' % scaf_idx
			
			# update index.
			idx += 1
			
		# update scaffold.
		scaf_idx += 1

	# return edges.
	return ne
			
########### script ################## 

# make graphs.
RG, ref_edges = make_agp(ref_file)
TG, test_edges = make_agp(test_file)

# get list of nodes.
nodes = [n for n in RG.nodes()]

# ensure all nodes present in TG.
TG, test_edges = ensure_consistent(TG, nodes, ref_edges, test_edges)

# build true adj set.
true_adj = build_adj(RG)

# orient the test graph components.
TG = orient_TG(nodes, true_adj, TG)

# calculate 4 parameters.
#A, B, C, D = compute_whole(nodes, true_adj, TG)
A, B, C, D, A_SET = compute_formula(nodes, true_adj, RG, TG)

# sanity check the 4 parameters.
sanity_check(RG, TG, A, B, C, D)

# calculate orientation.
orien_good, orien_bad = compute_orientation(ref_edges, test_edges, TG, true_adj)

# get gap deviation.
gap_dev = gap_deviation(A_SET, ref_edges, test_edges)

# get the runtime.
runtime = get_runtime(test_file)

# get the N50.
n50 = calculate_n50(test_edges)

# get the TP50.
ne = trueify_scaffold(test_edges, TG, A_SET)
tp50 = calculate_n50(ne)

# print the result.
res = []
res.append(name)
res.append(method)
res.append(orien_good)
res.append(orien_bad)
res.append(A)
res.append(B)
res.append(C)
res.append(D)
res.append(gap_dev)
res.append(runtime)
res.append(float(A) / float(A+C))
res.append(float(A) / float(A+B))
res.append(float((A*D)-(B*C)) / float(math.sqrt((A+B)*(A+C)*(D+B)*(D+C))))
res.append(n50)
res.append(tp50)


print '\t'.join([str(x) for x in res])
