'''
solves scaffolding problem using ILP
on pairs only. Then weighted bi-partite matching
to find the paths.
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
output_agp_file = sys.argv[2]
result_file = sys.argv[3]
cutoff = int(sys.argv[4])
previous_agp = sys.argv[5]
previous_res = sys.argv[6]

# hardcoded formatting.
input_dir = "%s/input" % base_dir
input_nodes_file = "%s/nodes.txt" % input_dir
input_edges_file = "%s/edges.txt" % input_dir
input_bundles_file = "%s/bundles.txt" % input_dir

graph_dir = "%s/graphs" % base_dir
decomp_0_file = "%s/decomp_zero.txt" % graph_dir
decomp_1_file = "%s/decomp_one.txt" % graph_dir
decomp_2_file = "%s/decomp_two.txt" % graph_dir

agp_dir = "%s/agps" % base_dir
cplex_log_file = "%s/logs/cplex.log" % base_dir
cplex_err_file = "%s/logs/cplex.err" % base_dir

int_dir = "%s/intermediate" % base_dir

#################### classes ####################

class GenState(object):
	bound = False

class SpqrSolve(object):
	''' solves ILP based on decomposition files '''
	
	def __init__(self, nodes, bundles, DG, ilp_obj, sol_obj, cutoff, node_force, bundle_force):
		''' sets tracking variables and starts solving '''
		
		# save variables.
		self._nodes = nodes
		self._bundles = bundles
		self._DG = DG
		self._ilp_obj = ilp_obj
		self._sol_obj = sol_obj
		self._cutoff = cutoff
		self._maxnidx = self._nodes[self._nodes.size - 1]['node_idx'] + 1
		
		# tracking
		self._loaded = False
		self._solved = False
		
		# register any presolved parts.
		self._pre_node = node_force
		self._pre_bundle = bundle_force
				
		# start decomposition.
		self._zero()
		
	def _load(self, active, biconts=[]):
		''' creates objects for solving '''
		
		# build lists.
		nlist = build_nlist(self._nodes, active)
		blist = build_blist(self._bundles, active)
		tlist = build_tlist(self._bundles, active)
		
		# build paritial solution.
		partial = PartialSolution(nlist.size, blist.size)
		
		# load the constraints.
		self._ilp_obj.load(nlist, blist, tlist, partial)
		
		
		# check if we need to force any orientation.
		isec = active.intersection(self._pre_node)
		if len(isec) > 0:
			for i in isec:
				# force orientation to be same.
				self._ilp_obj.force_one(i, self._pre_node[i])
		
		# check if we need to force any path.
		for i in range(blist.size):
			
			# check.
			idxa = blist[i]['idxa']
			idxb = blist[i]['idxb']
			
			# skip non-intersection.
			if idxa not in isec or idxb not in isec: continue
			
			# force if a hit.
			key = -1
			if (idxa, idxb) in self._pre_bundle:
				if self._pre_bundle[(idxa, idxb)][0] == 1:
					key = (idxa, idxb)
				elif self._pre_bundle[(idxa, idxb)][1] == 1:
					key = (idxb, idxa)
			
			elif (idxb, idxa) in self._pre_bundle:
				if self._pre_bundle[(idxb, idxa)][0] == 1:				
					key = (idxb, idxa)
				elif self._pre_bundle[(idxb, idxa)][1] == 1:
					key = (idxa, idxb)
			
			if key != -1:
				self._ilp_obj.force_path(key[0], key[1])

		# save the partial solution.
		self._partial = partial
		self._loaded = True

		
	def _force(self, i, val):
		''' forces orientation constraint '''
		
		# first check if is pre_solves.
		if i in self._pre_node:
			return
			
		# otherwise constrain.
		self._ilp_obj.force_one(i, val)
		
	def _clear(self):
		''' clears out a partial solution '''
		
		# clear solution.
		del self._partial
		
		# clear out ILP.
		self._ilp_obj.clear()
		
		# clear the var.
		self._loaded = False
		self._solved = False
		
	def _solve(self):
		''' executes the loaded partial solution '''
		
		# sanity.
		if self._loaded == False:
			logging.error("ilp not loaded")
			sys.exit(1)
		
		# call the solve method.
		self._partial, self._objval = self._ilp_obj.solve("cheese.lp")
		
		# note we finished.
		self._solved = True
		
		# return solution.
		return self._partial, self._objval
		
	def _intersection(self, active):
		''' finds intersection between solved and active '''
		
		# sanity.
		if self._solved == False:
			logging.error("bad intersection logic")
			sys.exit(1)			
			
		# retrieve intersection.
		return self._sol_obj.get_intersection(active)
		
	def _apply(self, mode=None):
		''' applies the solution '''
		
		# sanity.
		if self._solved == False:
			logging.error("bad apply logic")
			sys.exit(1)
		
		# call the apply method.
		if mode == None:
			self._sol_obj.apply_partial(self._partial)
			
		elif mode == "bi":
			self._sol_obj.apply_partial(self._partial)

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
			if len(s0) < self._cutoff:
				zero_smalls.append(n0)
			else:
				zero_larges.append(n0)
				
		# total sizes.
		tot = len(zero_smalls) + len(zero_larges)
				
		# solve the large components.
		for n0 in zero_larges:
			
			# simplify.
			g0 = self._DG.node[n0]['graph']
			s0 = self._DG.node[n0]['set']
			
			# solve the bucket.
			sol, val = self._one(g0, s0)
			
			# apply the solution.
			self._sol_obj.apply_partial(sol)
		
				
		# build the small buckets.
		buckets = []
		bucket = set()
		for n0 in zero_smalls:
			
			# simplify.
			s0 = self._DG.node[n0]['set']
			
			# add set to bucket.
			bucket = bucket.union(s0)
			
			# check if bucket is large enough.
			if len(bucket) > self._cutoff:
				
				#  add bucket to list.
				buckets.append(bucket)
				
				# clear the bucket.
				bucket = set()
		
		# add last one.
		buckets.append(bucket)
		
		# solve the small components.
		for bucket in buckets:
			
			# solve the final bucket.
			logging.info("solving zero: %d" % len(bucket))
			self._load(bucket)
			self._solve()
			self._apply()
			self._clear()

			
	def _one(self, g, s0):
		''' solves one component '''
		
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

		# choose root.
		root = g.nodes()[0]

		# solve bottom up.
		logging.info("solving one (bottom): %d" % len(s0))
		for p, q, kids in self._dfs_gen(g, root, order=False):
			
			# grab active set.
			s1 = g.node[p]['set']
			g1 = g.node[p]['graph']
			large = len(s1) > self._cutoff
			
			# if its the root just solve.
			if q == -1:
				
				# just solve (hope its not a 2-component).
				if large == False:		
					# solve normal.
					self._load(s1)
					sol, val = self._solve()
					self._clear()
				else:
					# solve 2 decomp.
					print "NOT WRITTEN"
					sys.exit()
				
				# save to root.
				g.node[p]['sol'] = sol
				g.node[p]['val'] = val
			
			else:
			
				# get cut.
				cut = g[p][q]['cut']
				
				# solve 2 combinations of orientations.
				logging.info("solving one: %d" % len(s1))
				for x in [0, 1]:
					
					# switch on solve.
					if large == False:
					
						# solve normal.
						self._load(s1)
						self._force(cut, x)
						sol, val = self._solve()
						self._clear()
						
					else:
												
						# solve decompose.
						sol, val = self._two(g1, s1, cut, x)
						
					# save info.
					g.node[p]['sol_%i' % x] = sol
					g.node[p]['val_%i' % x] = val
				
		# prepare a partial solution.
		nlist = build_nlist(self._nodes, s0)
		blist = build_blist(self._bundles, s0)
		tlist = build_tlist(self._bundles, s0)
		partial = PartialSolution(nlist.size, blist.size)
		val = 0.0
		
		# apply solution top down.
		logging.info("applying one (top): %d" % len(s0))
		for p, q, kids in self._dfs_gen(g, root, order=True):
			
			# grab active set.
			s1 = g.node[p]['set']
			
			# just apply root.
			if q == -1:
				partial.merge(g.node[p]['sol'])
				val += g.node[p]['val']
				continue
			
				
			# determine what the parent did.
			cut = g[p][q]['cut']
			o = partial.get_orien(cut)
				
			# apply that solution.
			partial.merge(g.node[p]['sol_%i' % o])
			val += g.node[p]['val_%i' % o]
			
		# return solution.
		return partial, val
					
			
	def _two(self, g, s1, cut1, val1):
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
		
		# solve bottom up.
		logging.info("solving two (bottom): %d" % len(s1))
		for p, q, kids in self._dfs_gen(g, root, order=False):
			
			# grab active set.
			s2 = g.node[p]['set']
			
			# if its the root just solve.
			if q == -1:
				
				# just solve.
				self._load(s2)
				if cut1 in s2:
					self._force(cut1, val1)
				sol, val = self._solve()
				self._clear()
				
				# save to root.
				g.node[p]['sol'] = sol
				g.node[p]['val'] = val
			
			else:
			
				# get cuts.
				cuts = g[p][q]['cut']
				
				# solve 4 combinations of orientations.
				for x, y in itertools.product([0,1], repeat=2):
					
					# solve.
					self._load(s2)
					if cut1 in s2:
						self._force(cut1, val1)
					self._force(cuts[0], x)
					self._force(cuts[1], y)
					sol, val = self._solve()
					self._clear()
					
					# save info.
					g.node[p]['sol_%i_%i' % (x,y)] = sol
					g.node[p]['val_%i_%i' % (x,y)] = val
				
		# prepare a partial solution.
		nlist = build_nlist(self._nodes, s1)
		blist = build_blist(self._bundles, s1)
		tlist = build_tlist(self._bundles, s1)
		partial = PartialSolution(nlist.size, blist.size)
		val = 0.0
							
		# apply solution top down.
		logging.info("applying two (top): %d" % len(s1))
		for p, q, kids in self._dfs_gen(g, root, order=True):
			
			# grab active set.
			s2 = g.node[p]['set']
			
			# just apply root.
			if q == -1:
				partial.merge(g.node[p]['sol'])
				val += g.node[p]['val']
				continue
			
				
			# determine what the parent did.
			cuts = g[p][q]['cut']
			oa = partial.get_orien(cuts[0])
			ob = partial.get_orien(cuts[1])
				
			# apply that solution.
			partial.merge(g.node[p]['sol_%i_%i' % (oa, ob)])
			val += g.node[p]['val_%i_%i' % (oa, ob)]
		
		# return the solution.
		return partial, val
		
		

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
		
def call_agp_gaps(agp_file, nodes):
	''' calls agp gaps'''
	
	# create node lookup.
	lookup = create_lookup(nodes)
	
	# load the agp array.
	agp_edges = load_agps(agp_file)
	
	# ensure sorted by scaffname and scafidx.
	agp_edges.sort(order=['scaf_name','scaf_idx'])
	
	# build list of component offsets.
	offsets = dict()
	for i in range(agp_edges.size):

		# skip non contigs.
		if agp_edges[i]['comp_type'] != "W": continue
		
		# record index.
		if agp_edges[i]['scaf_name'] not in offsets:
			offsets[agp_edges[i]['scaf_name']] = list()
		offsets[agp_edges[i]['scaf_name']].append(i)
		
		# add bundle info to this.
		gaps = dict()
		for key in offsets:
			
			# loop over edges.
			for i in range(len(offsets[key]) - 1):
				
				# get AGP edge.
				ea = agp_edges[offsets[key][i]]
				eb = agp_edges[offsets[key][i+1]]
				
				# get index.
				idxa = lookup[ea['comp_name']]
				idxb = lookup[eb['comp_name']]
				
				# get gap.
				gaps[(idxa,idxb)] = eb['scaf_start'] - ea['scaf_stop']
				
	return gaps
		
def load_previous(agp_file, nodes):
	''' loads info from previous AGP'''
	
	# create node lookup.
	lookup = create_lookup(nodes)
	
	# load the agp array.
	agp_edges = load_agps(agp_file)
	
	# ensure sorted by scaffname and scafidx.
	agp_edges.sort(order=['scaf_name','scaf_idx'])
	
	# build list of component offsets.
	orien = dict()
	offsets = dict()
	for i in range(agp_edges.size):

		# skip non contigs.
		if agp_edges[i]['comp_type'] != "W": continue
		
		# save orientation.
		orien[lookup[agp_edges[i]['comp_name']]] = agp_edges[i]['comp_orien']
		
		# record index.
		if agp_edges[i]['scaf_name'] not in offsets:
			offsets[agp_edges[i]['scaf_name']] = list()
		offsets[agp_edges[i]['scaf_name']].append(i)
		
	# add bundle info to this.
	gaps = dict()
	active = set()
	for key in offsets:
		
		# loop over edges.
		for i in range(len(offsets[key]) - 1):
			
			# get AGP edge.
			ea = agp_edges[offsets[key][i]]
			eb = agp_edges[offsets[key][i+1]]
			
			# get index.
			idxa = lookup[ea['comp_name']]
			idxb = lookup[eb['comp_name']]
			
			# get gap.
			gaps[(idxa,idxb)] = eb['scaf_start'] - ea['scaf_stop']
			
			# note its active.
			active.add((idxa,idxb))
			
	# return gaps, active set and node set.
	return gaps
		
########### script ################## 

# make sure output dir is there.
create_dir(agp_dir)
create_dir(int_dir)

# load hdf5 information.
logging.info("loading data arrays")
nodes = load_nodes(input_nodes_file)
edges = load_edges(input_edges_file)
bundles = load_bundles(input_bundles_file)

# identify AGP edges.
agp_gaps = {}
if previous_agp != "-":
	logging.info("loading agp edge distances")
	agp_gaps = load_previous(previous_agp, nodes)

# identify previous solution.
node_force = {}
bundle_force = {}
if previous_res != "-":
	logging.info("loading previous solution")
	fin = open(previous_res, "rb")
	for line in fin:
		# tokenize.
		tmp = line.strip().split()
		idxa = int(tmp[0])
		idxb = int(tmp[1])
		oriena = int(tmp[2])
		orienb = int(tmp[3])
		patha = int(tmp[4])
		pathb = int(tmp[5])
		
		# fill node force.
		node_force[idxa] = oriena
		node_force[idxb] = orienb
		
		# fill bundle.
		bundle_force[(idxa,idxb)] = (patha, pathb)
	fin.close()
		

# create solution object.
sol_obj = SpqrSolution(nodes, bundles)

logging.info("adding gap estimates")
gap_ests = gap_avg(nodes, edges, bundles, agp_gaps)

# create ILP object.
ilp_obj = SpqrIlp(cplex_log_file, cplex_err_file)

# load decomposition.
DG = load_decomposition(decomp_0_file, decomp_1_file, decomp_2_file)

# create an SPQR solve object.
solver = SpqrSolve(nodes, bundles, DG, ilp_obj, sol_obj, cutoff, node_force, bundle_force)

# finalize solution.
sol_obj.finalize()

# create a directed graph from solution.
TG = Directed(nodes, edges, bundles)
TG.create_edges(sol_obj)

# create a linear graph.
LG = Linear(nodes, edges)
LG.set_orientation(sol_obj)

# order potential edges.
LG.order_bundles(bundles, sol_obj.get_bundles())

# loop over eachs subgraph.
logging.info("solving for path")
for subg in TG.comps():
	
	# create the bipartite flow graph.
	FG = Flow(subg)
	
	# calculate path.
	elist = FG.pathify()
		
	# add to linear graph.
	LG.create_elist(elist)

# apply them.
LG.set_gaps(gap_ests)

# remove cycles.
logging.info("removing cycles")
LG.remove_cycles()

# end of execution.
time_stop = time.time()

# verify it.
LG.verify()

# write agp.
LG.write_agp(output_agp_file)

# append runtime to agp.
fin = open(output_agp_file, "rb")
lines = fin.readlines()
fin.close()

lines.append("RUNTIME:\t%f\n" % (time_stop - time_start))

fout = open(output_agp_file, "wb")
fout.write(''.join(lines))
fout.close()

logging.info("COMPLETED!")

# write out solution bundles to a file.
sb = sol_obj.get_bundles()
sn = sol_obj.get_nodes()

active = set(LG.edges())
fout = open(result_file, "wb")
for s in sb:
	if (s['idxa'], s['idxb']) in active:
		fout.write("%i\t%i\t%i\t%i\t%i\t%i\n" % (s['idxa'], s['idxb'], sn[s['idxa']]['orien'], sn[s['idxb']]['orien'], s['X'], s['Y']))
	elif (s['idxb'], s['idxa']) in active:
		fout.write("%i\t%i\t%i\t%i\t%i\t%i\n" % (s['idxa'], s['idxb'], sn[s['idxa']]['orien'], sn[s['idxb']]['orien'], s['X'], s['Y']))
fout.close()

	
	

