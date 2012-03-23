'''
verifies the decomposition graphs are valid.
'''
# program imports.
from data_structs.nodes import load_nodes
from data_structs.nodes import create_lookup
from data_structs.edges import load_edges
from data_structs.edges import save_edges
from data_structs.weights import load_weights
from data_structs.bundles import load_bundles
from data_structs.agps import load_agps
from data_structs import types

from optimize.ilp_pair_cplex import SpqrIlp
from optimize.ilp_pair_cplex import BiConstraints
from optimize.ilp_pair_cplex import nlist_dt
from optimize.ilp_pair_cplex import blist_dt
from optimize.ilp_pair_cplex import tlist_dt
from optimize.solution_pair import SpqrSolution
from optimize.solution_pair import PartialSolution

from graphs.Directed import Directed
from graphs.Flow import Flow
from graphs.Linear import Linear
from graphs.BiConnected import bicomponents

from gap_estimate.average import gap_avg

# system imports.
from collections import deque
import random
import time
import numpy as np
import subprocess
import itertools
import networkx as nx
import logging
import sys
import os

time_start = time.time()

# logging.
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )
#logging.basicConfig(level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )

# parameters.
base_dir = os.path.abspath(sys.argv[1])

# hardcoded formatting.
input_dir = "%s/input" % base_dir
input_nodes_file = "%s/nodes.txt" % input_dir
input_edges_file = "%s/edges.txt" % input_dir
input_bundles_file = "%s/bundles.txt" % input_dir

graph_dir = "%s/graphs" % base_dir
decomp_0_file = "%s/decomp_zero.txt" % graph_dir
decomp_1_file = "%s/decomp_one.txt" % graph_dir
decomp_2_file = "%s/decomp_two.txt" % graph_dir


#################### classes ####################

class TestSol(object):
	''' stores flag if variable was solved '''
	
	def __init__(self):
		self.visited = set()
	
	def mark(self, sol_1, sol_2, size):
		''' saves partial solution '''

		# check for intersection.
		z = sol_1.intersection(sol_2)
		
		# check for intersection.
		if len(z) != size:
			logging.error("BAD INTERSECTION")
			print sol_1, sol_2, z, size
			for x in self.visited:
				print x
			
			sys.exit(1)

		# make the solution.
		self.visited = self.visited.union(sol)
		
		

class SpqrTest(object):
	''' solves ILP based on decomposition files '''
	
	def __init__(self, nodes, bundles, DG, sol_obj):
		''' sets tracking variables and starts solving '''
		
		# save variables.
		self._nodes = nodes
		self._bundles = bundles
		self._DG = DG
		self._sol_obj = sol_obj
				
		# start decomposition.
		self._zero(DG)

	def _zero(self, G):
		''' bootstraps the solution process ''' 

		# loop over large nodes.
		for n in G.nodes():
			
			# skip n.
			if G.node[n]['graph'] == None:
				continue
			
			# loop into this.
			self._solve(G.node[n]['graph'])
			
	def _solve(self, G):
		''' solves the decomposition '''
		
		#'''
		solved = set()
		for n in G.nodes():
			
			# check for intersection.
			inter = solved.intersection(G.node[n]['set'])
			if len(inter) > 1:
				print n, inter
				print "WAT"
				sys.exit()
			
			solved = solved.union(G.node[n]['set'])
		#'''	
		return	
		
		# bottom up solution.
		for p, parent, isleaf, isroot in G.graph['bottom']:
			
			# get set.
			s0 = G.node[p]['set']
			s1 = G.node[parent]['set']
			
			# check for intersection.
			cut = G[p][parent]['cut']
			
			# find intersection size.
			cut_size = len(cut)
			if isleaf == True:
				cut_size = 0
			
			# solve it.
			self._sol_obj.mark(s0, s1, cut_size)

		

#################### functions ####################

def draw(name, subg):
	''' draw an interactive graph.'''
	# write dot file.
	dot = "./tmp.dot"
	pic = "/home/jlindsay/central/transfer/%s.jpg" % name
	
	nx.write_dot(subg, dot)
	
	# convert to picture using neato.
	subprocess.call(["neato", "-Tjpg", "-o", pic, dot])
	
	# remove dot.
	subprocess.call(["rm","-f",dot])

def create_dir(file_path):
	''' creates a dir. '''
	if os.path.isdir(file_path) == False:
		subprocess.call(["mkdir",file_path])

def load_decomposition2(decomp_0_file, decomp_1_file, decomp_2_file):
	''' loads skeleton from decomposition file.'''
	
	# build dictionary of decomposition sets.
	DG = nx.Graph()
	
	# loop over all files.
	#for file_path in [decomp_0_file, decomp_1_file, decomp_2_file]:
	for file_path in [decomp_0_file, decomp_1_file]:
		
		# loop over all lines.
		fin = open(file_path, "rb")
		for line in fin:
			
			# tokenize.
			line = line.strip().strip(",").split("\t")
			stage_type = int(line[0].replace("stage=",""))
			line_type = line[1].replace("type=","")

			# look for nodes.
			if line_type == "N":
			
				# get index.
				tmp = line[2].replace("idx=","").strip().strip(",").split(",")
				idx0 = int(tmp[0])
				idx1 = -1
				idx2 = -1
				
				if idx1 > 25: continue
				
				if len(tmp) > 1:
					idx1 = int(tmp[1])
				if len(tmp) > 2:
					idx2 = int(tmp[2])

				# build list of active.
				tmp1 = line[3].replace("set=","").strip(",")
				tmp2 = tmp1.split(",")			
				
				# add this set to list.
				tmp = set()
				for x in tmp2:
					tmp.add(int(x))

				# add set to graph.
				if idx0 != -1 and  DG.has_node(idx0) == False:
					DG.add_node(idx0, {'graph':nx.Graph(), 'set':set()})

				if idx1 != -1 and DG.node[idx0]['graph'].has_node(idx1) == False:
					DG.node[idx0]['graph'].add_node(idx1, {'graph':nx.Graph(), 'set':set()})

				if idx2 != -1 and DG.node[idx0]['graph'].node[idx1]['graph'].has_node(idx2) == False:
					DG.node[idx0]['graph'].node[idx1]['graph'].add_node(idx2, {'graph':None, 'set':set()})
				
				# add set at maximum key size.
				if idx0 != -1 and idx1 != -1 and idx2 != -1:
					DG.node[idx0]['graph'].node[idx1]['graph'].node[idx2]['set'] = tmp
					
				elif idx0 != -1 and idx1 != -1:
					DG.node[idx0]['graph'].node[idx1]['set'] = tmp
					
				elif idx0 != -1:
					DG.node[idx0]['set'] = tmp

		fin.close()
	
	# clear out non-decomposed entries.
	recurse_clear(DG)
	
	# compute edges.
	recurse_edgeify(DG)
	
	# compute decomposition paths.
	for n in DG.nodes():
		
		# blowout this node if necessary.
		if DG.node[n]['graph'] != None:
			blowout(DG.node[n]['graph'])
			
			
		G = DG.node[n]['graph']
		
		
		solved = set()
		for n in G.nodes():
			
			# check for intersection.
			inter = solved.intersection(G.node[n]['set'])
			if len(inter) > 1:
				print n, inter
				print "WAT"
				sys.exit()
			
			solved = solved.union(G.node[n]['set'])
	
	sys.exit()
	
	# return the final object.
	return DG
	
def blowout(G):
	''' builds ordered edge list. '''
	
	# bootstrap edges.
	root = G.nodes()[0]
	Vnew = set([root])
	V = set(G.nodes())
	
	# ensure its connected.
	if nx.is_connected(G) == False:
		logging.error("shouldn't be here")
		sys.exit(1)
	
	# grow the sets.
	elist = list([root, -1])
	while V != Vnew:
		
		# loop over adjacencies.
		nlist = list(Vnew)
		for n1 in nlist:
			for n2 in G.neighbors(n1):
				if n2 in Vnew: continue
				
				# further blowout this node.
				if G.node[n2]['graph'] != None:
					blowout(G.node[n2]['graph'])
				
				# add edge.
				elist.append( (n2, n1) )
				
				# add node to set.
				Vnew.add(n2)

	# encode this.
	G.graph['elist'] = elist
	
	
def dfs_bottom(G, tr=False):
	''' yields [node, parent, isleaf(bool), isroot(bool)] bottom up manner '''
	
	# bootstrap DFS.
	if tr == False:
		root = G.nodes()[0]
	else:
		root = tr
	stack = list()
	parent = dict()
	children = dict()
	visited = set()
		
	stack.append(root)
	parent[root] = - 1
	parent = {root:-1}
	
	# begin DFS.
	while len(stack) > 0:
		
		# peek at end.
		p = stack[-1]
				
		# mark as visited.
		visited.add(p)
		
		# check children.
		cnt = 0
		for q in G.neighbors(p):
			
			# skip visited.
			if q in visited: continue
			
			# add to stack.
			stack.append(q)
			cnt += 1
			
			# note parent.
			parent[q] = p
			
			# note children.
			if p not in children:
				children[p] = list()
			children[p].append(q)
		
		# continue if all kids not processed.
		if cnt != 0: 
			continue

		# remove from stack.
		p = stack.pop()
		
		# get leaf status.
		if p in children:
			isleaf = False
		else:
			isleaf = True
			
		# get root status.
		if parent[p] == -1:
			isroot = True
		else:
			isroot = False
		
		# yield the info.
		yield p, parent[p], isleaf, isroot
		
def dfs_top(G):
	''' yields [node, parent, isleaf(bool), isroot(bool)] top down manner '''
	
	# bootstrap DFS.
	root = G.nodes()[0]
	stack = list()
	parent = dict()
	children = dict()
	visited = set()
		
	stack.append(root)
	parent[root] = - 1
	parent = {root:-1}
	
	# begin DFS.
	while len(stack) > 0:
		
		# peek at end.
		p = stack[-1]
				
		# mark as visited.
		visited.add(p)
		
		# check children.
		cnt = 0
		for q in G.neighbors(p):
			
			# skip visited.
			if q in visited: continue
			
			# add to stack.
			stack.append(q)
			cnt += 1
			
			# note parent.
			parent[q] = p
			
			# note children.
			if p not in children:
				children[p] = list()
			children[p].append(q)
		
		# check if this node has parents.
		if parent[p] == -1:
			isroot = True
		else:
			isroot = False
			
		# check if this node is a leaf.
		if p in children:
			isleaf = False
		else:
			isleaf = True
			
		# yield this was encountered.
		yield p, parent[p], isleaf, isroot
		
		# continue if all kids not processed.
		if cnt != 0: 
			continue

		# remove from stack.
		p = stack.pop()
		
	
def recurse_paths(G):
	''' replaces each graph with a directed version '''
	
	# create nlist.
	nlist = G.nodes()
		
	# digg deeper for each node.
	for n in nlist:
			
		# check if we add graph.
		if G.node[n]['graph'] != None:
			recurse_paths(G.node[n]['graph'])
					
					
	# build the edge lists.
	G.graph['bottom'] = list()
	for x in dfs_bottom(G):
		if x[1] != -1:
			G.graph['bottom'].append(x)

	G.graph['top'] = list()
	for x in dfs_top(G):
		if x[1] != -1:
			G.graph['top'].append(x)
		

	
def recurse_edgeify(G):
	''' adds edges to graph based on overlap '''
	
	# create nlist.
	nlist = G.nodes()
	
	# loop over nodes of graph.
	for i in range(len(nlist)):
		
		# check if this requires further edgeification.
		if G.node[nlist[i]]['graph'] != None:
			recurse_edgeify(G.node[nlist[i]]['graph'])
			
		# look for edges.
		for j in range(i+1, len(nlist)):
			
			# check for intersection.
			inter = G.node[nlist[i]]['set'].intersection(G.node[nlist[j]]['set'])
			
			if len(inter) > 2:
				print "SAAAD"
			
			# if intersection annotate it and add edge.
			if len(inter) > 0:
				G.add_edge(nlist[i], nlist[j], {'cut': inter})
				
	
def load_decomposition(decomp_0_file, decomp_1_file, decomp_2_file):
	''' loads skeleton from decomposition file.'''
	
	# read zero file.
	fin = open(decomp_0_file, "rb")
	zero_lines = [x.strip().split("\t") for x in fin.readlines()]
	fin.close()
	
	# read one file.
	fin = open(decomp_1_file, "rb")
	one_lines = [x.strip().split("\t") for x in fin.readlines()]
	fin.close()	
	
	# read two file.
	fin = open(decomp_2_file, "rb")
	two_lines = [x.strip().split("\t") for x in fin.readlines()]
	fin.close()	

	# create master graph.
	G = nx.DiGraph()
	
	##### zero #####
	
	# add zero skeleton nodes.
	for line in zero_lines:
		
		# sanity check.
		if line[0] != "stage=0":
			logging.error("stage 0 error: 1")
			sys.exit(1)

		if line[1].count("type=") == 0:
			logging.error("stage 0 error: 2")
			sys.exit(1)

		# tokenize.
		line_type = line[1].replace("type=","")

		# switch on type.
		if line_type == "N":
		
			# get index.
			tmp = line[2].replace("idx=","").split(",")
			node_idx = int(tmp[0])
			
			# build list of active.
			tmp1 = line[3].replace("set=","").strip(",")
			tmp2 = tmp1.split(",")
			
			nset = set()
			for x in tmp2:
				nset.add(int(x))
			
			# add node.
			G.add_node( node_idx, {\
				'set': nset,\
				'graph': nx.Graph()\
			})
		
	##### one #####
	
	# add one skeleton nodes.
	for line in one_lines:
		
		# sanity check.
		if line[0] != "stage=1":
			logging.error("stage 1 error: 1")
			sys.exit(1)

		if line[1].count("type=") == 0:
			logging.error("stage 1 error: 2")
			sys.exit(1)

		# tokenize.
		line_type = line[1].replace("type=","")

		# switch on type.
		if line_type == "N":
		
			# get index.
			tmp = line[2].replace("idx=","").split(",")
			idx0 = int(tmp[0])
			idx1 = int(tmp[1])
			
			# build list of active.
			tmp1 = line[3].replace("set=","").strip(",")
			tmp2 = tmp1.split(",")
			
			nset = set()
			for x in tmp2:
				nset.add(int(x))
			
			# add node.
			G.node[idx0]['graph'].add_node( idx1, {\
				'set': nset,\
				'graph': nx.Graph()\
			})
		elif line_type == "E":

			# get index.
			tmp = line[2].replace("idx1=","").split(",")
			idxa0 = int(tmp[0])
			idxa1 = int(tmp[1])
			
			tmp = line[3].replace("idx2=","").split(",")
			idxb0 = int(tmp[0])
			idxb1 = int(tmp[1])
			
			# cut cut.
			cut = int(line[4].replace("cut=",""))
			
			# add edges.
			G.node[idx0]['graph'].add_edge( idxa1, idxb1, {\
				'cut': cut
			})
		elif line_type == "R":
			
			# get index.
			tmp = line[2].replace("idx=","").split(",")
			idx0 = int(tmp[0])
			root = int(tmp[1])
			
			# add root to tree.
			G.node[idx0]['graph'].graph['root'] = root
			
	##### two #####
	
	# add two skeleton nodes.
	for line in two_lines:
		
		# sanity check.
		if line[0] != "stage=2":
			logging.error("stage 2 error: 1")
			sys.exit(1)

		if line[1].count("type=") == 0:
			logging.error("stage 2 error: 2")
			sys.exit(1)

		# tokenize.
		line_type = line[1].replace("type=","")

		# switch on type.
		if line_type == "N":
		
			# get index.
			tmp = line[2].replace("idx=","").split(",")
			idx0 = int(tmp[0])
			idx1 = int(tmp[1])
			idx2 = int(tmp[2])
			
			# build list of active.
			tmp1 = line[3].replace("set=","").strip(",")
			tmp2 = tmp1.split(",")
			
			nset = set()
			for x in tmp2:
				nset.add(int(x))
			
			# add node.
			G.node[idx0]['graph'].node[idx1]['graph'].add_node( idx2, {\
				'set': nset,\
				'graph': nx.Graph()\
			})
		elif line_type == "E":

			# get index.
			tmp = line[2].replace("idx1=","").split(",")
			idxa0 = int(tmp[0])
			idxa1 = int(tmp[1])
			idxa2 = int(tmp[2])
			
			tmp = line[3].replace("idx2=","").split(",")
			idxb0 = int(tmp[0])
			idxb1 = int(tmp[1])
			idxb2 = int(tmp[2])
			
			# sanity.
			if idxa0 != idxb0 or idxa1 != idxb1:
				logging.error("bad mojo")
				sys.exit()
				
			# simplify.
			idx0 = idxa0
			idx1 = idxa1
			
			# cut cut.
			tmp = line[4].replace("cut=","").strip(",").split(",")
			cut1 = int(tmp[0])
			cut = (cut1)
			
			if len(tmp) > 1:
				cut2 = int(tmp[1])
				cut = (cut1, cut2)
			
			# add edges.
			G.node[idx0]['graph'].node[idx1]['graph'].add_edge( idxa2, idxb2, {\
				'cut': cut
			})
		elif line_type == "R":
			
			# get index.
			tmp = line[2].replace("idx=","").split(",")
			idx0 = int(tmp[0])
			idx1 = int(tmp[1])
			root = int(tmp[2])
			
			# add root to tree.
			G.node[idx0]['graph'].node[idx1]['graph'].graph['root'] = root
		
	# return the graph.
	return G
		
def recurse_clear(G):
	for n in G.nodes():
		if G.node[n]['graph'] != None:
			if len(G.node[n]['graph'].nodes()) == 0:
				G.node[n]['graph'] = None
			else:
				recurse_clear(G.node[n]['graph'])
		
########### script ################## 

# load hdf5 information.
logging.info("loading data arrays")
nodes = load_nodes(input_nodes_file)
bundles = load_bundles(input_bundles_file)

# select a small example.
G = nx.Graph()
for x in bundles:
	G.add_edge(x['ctg_a_idx'], x['ctg_b_idx'])
	
for comp in nx.connected_components(G):
	if len(comp) > 25 and len(comp) < 100:
		active = comp
		break

G = G.subgraph(active)
draw("G", G)

# get bicomps and cuts.
bicomps, cuts = bicomponents(G)
bicomps = list(bicomps)

# map node to component.
nmap = dict()
for i in range(len(bicomps)):
	for n in bicomps[i]:
		if n not in nmap:
			nmap[n] = list()
		nmap[n].append(i)

# build BC graph.
BC = nx.Graph()

# add all cuts and bicomps to the graph.
idx = 0
cnodes = list()
bnodes = list()
for c in cuts:
	BC.add_node(idx, {'type':'C', 'cut':c})
	cnodes.append(idx)
	idx += 1
for bicomp in bicomps:
	BC.add_node(idx, {'type':'B', 'set':bicomp})
	bnodes.append(idx)
	idx += 1

# add edges from cuts to components. 
for c in cnodes:
	for b in bnodes:
		if BC.node[c]['cut'] in BC.node[b]['set']:
			BC.add_edge(c,b)
draw("BC", BC)


'''
# get edges from random cut root.
root = cnodes[0]
edges = [x for x in nx.dfs_edges(BC, root)]

# build tree.
T = nx.DiGraph()
for e in edges:
	T.add_edge(e[0], e[1])
draw("T", T)

# reverse the edges so we have bottom up.
edges.reverse()

for n in nx.dfs_postorder_nodes(BC):
	print n
sys.exit()
'''

# compute orderings.
root = cnodes[0]
preds = nx.dfs_predecessors(BC, root)
nlist = [x for x in nx.dfs_postorder_nodes(BC, root)]

# test solving.
solved = set()
for p in nlist:
	
	# get parent.
	parent = preds[p]
	
	# flip if is cut node.
	t = BC.node[p]['type']
	if t == 'C':
		continue
		
	# simplify.
	s0 = BC.node[p]['set']
	
	# check for intersection.
	inter = s0.intersection(solved)
	
	print p, parent, s0, s1
	
	# size should be 1 or 0.
	if len(inter) > 1:
		print "why so big, sad panda"
		print inter
		print solved
		sys.exit(1)
			
	# solve this.
	solved = solved.union(s0)
	


sys.exit()
A = nx.DiGraph()

for e in edges:
	A.add_edge(e[0], e[1])

# test




sys.exit()

solved = set()
for n in bicomps:
	
	# check for intersection.
	inter = solved.intersection(n)
	if len(inter) > 1:
		print n, inter, solved
		print "WAT"
		sys.exit()
	
	solved = solved.union(n)
	
	


sys.exit()


# test solution.
sol_test = TestSol()

# load decomposition.
logging.info("loading decomposition")
DG = load_decomposition2(decomp_0_file, decomp_1_file, decomp_2_file)

# create an SPQR solve object.
logging.info("solve the object")
solver = SpqrTest(nodes, bundles, DG, sol_test)


