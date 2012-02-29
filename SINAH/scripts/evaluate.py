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
	
	# build sorted edge set.
	edges = set()
	for i in range(agp_edges.size):
		
		# skip contigs themselves.
		if agp_edges[i]['comp_type'] != 'N': continue
		
		# add sorted edges.
		key = tuple(sorted([agp_edges[i-1]['comp_name'], agp_edges[i+1]['comp_name']]))
		edges.add(key)
		
	
	return edges, agp_edges

	
def compute_formula(true_adj, test_adj, n):
	''' computes statistics when TG is already oriented. This
	method computes A and B exactly, but uses formulas
	to find C and D.'''
	
	# set vars.
	m = len(true_adj)
	N = n * (n-1)

	# compute the A statistic set.
	A = true_adj.intersection(test_adj)
		
	# compute the B statistic set.
	B = test_adj.difference(true_adj)
	
	# compute the C statistic set.
	C = m - len(A)

	# compute the D statistic set.
	D = N - m - len(B)
		
	# return categories.
	return len(A), len(B), C, D, A

				
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
		if edges[i]['comp_type'] == 'W':
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
			
def trueify_scaffold(edges, a_set):
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
true_adj, ref_edges = make_agp(ref_file)
test_adj, test_edges = make_agp(test_file)

# counter number of contigs.
n = 0
for i in range(ref_edges.size):
	if ref_edges[i]['comp_type'] == 'W':
		n += 1

# calculate 4 parameters.
A, B, C, D, A_SET = compute_formula(true_adj, test_adj, n)

# get gap deviation.
gap_dev = gap_deviation(A_SET, ref_edges, test_edges)

# get the runtime.
runtime = get_runtime(test_file)

# get the N50.
n50 = calculate_n50(test_edges)

# get the TP50.
ne = trueify_scaffold(test_edges, A_SET)
tp50 = calculate_n50(ne)

# print the result.
res = []
res.append(name)
res.append(method)
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