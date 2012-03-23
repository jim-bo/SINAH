'''
merges small components into sizes
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
cutoff = int(sys.argv[2])

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
		

class OptimizeDecomp(object):
	''' runs through decomposition merging small components. '''
	
	def __init__(self, DG, cutoff):
		''' sets tracking variables and starts solving '''
		
		# save variables.
		self._DG = DG
		self._cutoff = cutoff
						
		# start decomposition.
		self._zero()
		
	def _dfs_gen(self, g, root, order=False):
		''' generator for DFS traversal '''
		
		# bootstrap a DFS.
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
			
			# yield in order.
			if p not in visited and order == True:
				if p in children:	kids = children[p]
				else:	kids = []
				yield p, parent[p], kids	
			
			# mark as visited.
			visited.add(p)
			
			# check children.
			cnt = 0
			for q in g.neighbors(p):
				
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

			# get children info.
			if p in children:
				kids = children[p]
			else:
				kids = []
			
			# continue if all kids not processed.
			if cnt != 0: 
					
				# continue
				continue
			
			# remove from stack.
			p = stack.pop()

			
			# yield at the bottom.
			if order == False:
				yield p, parent[p], kids				
		
	def _zero(self):
		''' solves large and small components '''

		# segregate components.
		zero_smalls = []
		zero_larges = []
		for n0 in self._DG.nodes():
			
			# simplify.
			g0 = self._DG.node[n0]['graph']
			s0 = self._DG.node[n0]['set']			
			
			# small / large check.
			if len(s0) == 0 or g0 == None:
				zero_smalls.append(n0)
			else:
				zero_larges.append(n0)
				
		# merge the small buckets.
		bucket = set()
		idx = 0
		for n0 in zero_smalls:
			
			# simplify.
			s0 = self._DG.node[n0]['set']
			
			# check if we add or finish.
			if len(s0) + len(bucket) > self._cutoff:
				
				# finalize new bucket.
				self._DG.add_node(idx, {'set':bucket, 'graph':None})
				
				# clear out bucket and idx.
				bucket = set()
				idx += 1
			
			# add set to bucket.
			bucket = bucket.union(s0)
			
			# remove old node.
			self._DG.remove_node(n0)
							
		# solve the large components.
		for n0 in zero_larges:
			
			# simplify.
			g0 = self._DG.node[n0]['graph']
			s0 = self._DG.node[n0]['set']
			
			if len(s0) == 0:
				sys.exit("say what")
				sys.exit(1)
			
			# merge components.
			self._one(g0, s0)

			
	def _one(self, g, s0):
		''' solves one component '''
		
		# sanity check.
		if len(s0) == 0 and len(g.nodes()) != 0:
			logging.error("bad sanity")
			sys.exit(1)
		
		# greedy graph constraction.
		logging.info("one contracting")
		psz = len(g.nodes())
		while 1 == 1:
			test = False
			for e in g.edges():
				
				# graph active sets.
				sa = g.node[e[0]]['set']
				sb = g.node[e[1]]['set']
				
				# check if we cant merge.
				if len(sa) + len(sb) > self._cutoff: continue
					
				# make new node info.
				n = "%s_%s" % (str(e[0]), str(e[1]))
				s = g.node[e[0]]['set'].union(g.node[e[1]]['set'])			
				
				# insert into graph.
				g.add_node(n, {'set':s, 'graph':None})
				
				# rehook edges.
				for q in g.neighbors(e[0]):
					if q == e[1]: continue
					cut = g[e[0]][q]['cut']
					g.remove_edge(e[0], q)
					g.add_edge(n, q, {'cut':cut})

				for q in g.neighbors(e[1]):
					if q == e[0]: continue
					cut = g[e[1]][q]['cut']
					g.remove_edge(e[1], q)
					g.add_edge(n, q, {'cut':cut})		
					
				# do contraction.
				g.remove_edge(e[0], e[1])
				g.remove_node(e[0])
				g.remove_node(e[1])
					
				# break and reset.
				test = True
				break
				
			# check if there was no mod.
			if test == False:
				break				
		logging.info("one contracted: %i to %i" % (psz, len(g.nodes())))

		# further simplify.
		for n in g.nodes():
			
			# simplify.
			s1 = g.node[n]['set']
			g1 = g.node[n]['graph']
			
			# check if has more.
			if g1 != None and len(g1.nodes()) > 0:
				
				# do two decomposition.
				self._two(g1, s1)
					
			
	def _two(self, g, s1):
		''' solves 2-decomposition '''
	
		# get root.
		root = g.graph['root']
	
		# greedy graph constraction.
		logging.info("two contracting")
		psz = len(g.nodes())
		while 1 == 1:
			test = False
			for e in g.edges():
				
				# graph active sets.
				sa = g.node[e[0]]['set']
				sb = g.node[e[1]]['set']
				
				# check if we cant merge.
				if len(sa) + len(sb) > self._cutoff: continue
				
				# dont touch root.
				if e[0] == root or e[1] == root: continue
					
				# make new node info.
				n = "%s_%s" % (str(e[0]), str(e[1]))
				s = g.node[e[0]]['set'].union(g.node[e[1]]['set'])			
				
				# insert into graph.
				g.add_node(n, {'set':s, 'graph':None})
				
				# rehook edges.
				for q in g.neighbors(e[0]):
					if q == e[1]: continue
					cut = g[e[0]][q]['cut']
					g.remove_edge(e[0], q)
					g.add_edge(n, q, {'cut':cut})

				for q in g.neighbors(e[1]):
					if q == e[0]: continue
					cut = g[e[1]][q]['cut']
					g.remove_edge(e[1], q)
					g.add_edge(n, q, {'cut':cut})		
					
				# do contraction.
				g.remove_edge(e[0], e[1])
				g.remove_node(e[0])
				g.remove_node(e[1])
					
				# break and reset.
				test = True
				break
				
			# check if there was no mod.
			if test == False:
				break				
		logging.info("two contracted: %i to %i" % (psz, len(g.nodes())))
		

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
		
def build_nlist(nodes, active):
	''' build nlist from active offset '''
		
	# allocate array.
	nlist = np.zeros(len(active), dtype=nlist_dt)
	
	# fill array.
	idx = 0
	for i in active:
		nlist[idx]['idx'] = nodes[i]['node_idx']
		idx += 1
	
	# return array.
	return nlist
	
def build_blist(bundles, active):
	''' build blist from active offset '''
	
	# count number of bundles.
	bsz = 0
	bidxs = []
	for i in range(bundles.size):
		if bundles[i]['ctg_a_idx'] in active and bundles[i]['ctg_b_idx'] in active:
			bsz += 1
			bidxs.append(i)
	
	# allocate array.
	blist = np.zeros(bsz, dtype=blist_dt)
	
	# fill array.
	idx = 0
	for i in bidxs:
		blist[idx]['idxa'] = bundles[i]['ctg_a_idx']
		blist[idx]['idxb'] = bundles[i]['ctg_b_idx']
		blist[idx]['WT_A'] = bundles[i]['WT_A']
		blist[idx]['WT_B'] = bundles[i]['WT_B']
		blist[idx]['WT_C'] = bundles[i]['WT_C']
		blist[idx]['WT_D'] = bundles[i]['WT_D']
		idx += 1

	# return array.
	return blist
	
def build_tlist(bundles, active):
	''' build blist from active offset '''
	
	# build dict of neighbors.
	neibs = dict()
	bidxs = []
	for i in range(bundles.size):
		
		# simplify.
		idxa = bundles[i]['ctg_a_idx']
		idxb = bundles[i]['ctg_b_idx']
		
		# skip.
		if idxa not in active or idxb not in active:
			continue
		
		# track it.
		bidxs.append(i)
		
		# bootstrap sets.
		if idxa not in neibs:
			neibs[idxa] = set()
		if idxb not in neibs:
			neibs[idxb] = set()
			
		# add to sets.
		neibs[idxa].add(idxb)
		neibs[idxb].add(idxa)
	
	# find triangles by edge.
	tris = set()
	for i in bidxs:
		
		# simplify.
		idxa = bundles[i]['ctg_a_idx']
		idxb = bundles[i]['ctg_b_idx']
		
		# build intersection.
		isec = neibs[idxa].intersection(neibs[idxb])
		
		# enumerate triangles.
		for idxc in isec:
			
			# make key.
			key = tuple( sorted( [idxa, idxb, idxc] ) )
			
			# insert to triangles.
			tris.add(key)
						
	# build array.
	tlist = np.zeros(len(tris), dtype=tlist_dt)
	
	# populate array.
	idx = 0
	for key in tris:
		tlist[idx]['idxa'] = key[0]
		tlist[idx]['idxb'] = key[1]
		tlist[idx]['idxc'] = key[2]
		idx += 1

	# return array.
	return tlist
		
def mark_used(edges, rnd):
	''' marks valid edges as used '''
	
	# loop over edges.
	for i in range(edges.size):
		
		# test valid.
		if edges[i]['invalid'] == 0:
			edges[i]['used'] = rnd
			
	# return edges.
	return edges
		
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

	# create neigbor lookup.
	

	# create master graph.
	G = nx.Graph()
	
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
			tmp = line[4].replace("cut=","").split(",")
			cut1 = int(tmp[0])
			cut2 = int(tmp[1])
			
			# add edges.
			G.node[idx0]['graph'].node[idx1]['graph'].add_edge( idxa2, idxb2, {\
				'cut': (cut1,cut2)
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
		
########### script ################## 

# load decomposition.
logging.info("loading decomposition")
DG = load_decomposition(decomp_0_file, decomp_1_file, decomp_2_file)

# create an SPQR solve object.
logging.info("solve the object")
o = OptimizeDecomp(DG, cutoff)


