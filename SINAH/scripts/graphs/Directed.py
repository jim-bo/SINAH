'''
builds a linear graph from solution.
'''
# program imports.
from data_structs.agps import agp_dt

# system imports.
from operator import itemgetter
import time
import numpy as np
import subprocess
import itertools
import networkx as nx
import logging
import sys
import os


############ class ###########
class Directed(nx.DiGraph):
	
	def __init__(self, nodes, edges, bundles):
		
		# call base class constructor.
		nx.DiGraph.__init__(self)
		
		# save pointers to important structures.
		self._nodes = nodes
		self._edges = edges
		self._bundles = bundles
		
		# add all the nodes right away.
		self._add_nodes()
	
	
	def comps(self):
		''' generator for connected components '''
		
		# loop over components.
		for comp in nx.weakly_connected_components(self):
			
			# create subgraph.
			subg = nx.DiGraph()
			
			# build node list.
			nlist = []
			for n in comp:
				nlist.append( (n, self.node[n]) )
				
			# add nodes.
			subg.add_nodes_from(nlist)
			
			# build edge list.
			elist = []
			for e in self.edges(comp):
				elist.append( (e[0], e[1], self[e[0]][e[1]]) )
				
			# add edges.
			subg.add_edges_from(elist)
			
			# yield the subgraph.
			yield subg
	
	def create_edges(self, sol, special=set()):
		''' creates graph from path variables '''
		
		# get nodes and bundles.
		sol_nodes = sol.get_nodes()
		sol_bundles = sol.get_bundles()
		
		# add the orientation to nodes.
		for i in range(sol_nodes.size):
			
			# simplify.
			idx = sol_nodes[i]['idx']
			orien = sol_nodes[i]['orien']

			# sanity check.
			if idx == -1:
				logging.error("missing orientation: %i" % idx)
				sys.exit(1)

			# set value.
			self.node[idx]['orien'] = orien
			
		
		# adds the bundles based on X,Y.
		to_add = []
		for i in range(sol_bundles.size):
			
			# get index.
			idxa = sol_bundles[i]['idxa']
			idxb = sol_bundles[i]['idxb']
			
			# sanity check non-special.
			if idxa not in special and idxb not in special:
				if idxa != self._bundles[i]['ctg_a_idx'] or idxb != self._bundles[i]['ctg_b_idx']:
					print sol_bundles[i]
					print self._bundles[i]
					logging.error("didn't think about that one did ya?")
					sys.exit(1)
			
			# skip empty.
			if sol_bundles[i]['X'] + sol_bundles[i]['Y'] == 0:
				continue
			
			# get weight of non special.
			if idxa not in special and idxb not in special:
				wt = 0.0
				for x in ['A', 'B', 'C', 'D']:
					if sol_bundles[i][x] == 1:
						wt = self._bundles[i]['WT_%s' % x]
				if wt == 0.0:
					continue
				wt = int(wt)
			else:
				wt = 9999999.99
			
			# add edge based on status.
			if sol_bundles[i]['X'] == 1:
				to_add.append( (idxa, idxb, {\
					'id':i,\
					'dist':None,\
					'weight':wt\
				}) )
			elif sol_bundles[i]['Y'] == 1:
				to_add.append( (idxb, idxa, {\
					'id':i,\
					'dist':None,\
					'weight':wt\
				}) )
		
		# add the edges at once.
		self.add_edges_from(to_add)
	
		
	def _add_nodes(self):
		''' adds all nodes no matter what '''
		
		# loop over nodes.
		to_add = []
		for i in range(self._nodes.size):
			
			# add node.
			to_add.append( (self._nodes[i]['node_idx'], {\
				'size': self._nodes[i]['ctg_width'], \
				'name': self._nodes[i]['ctg_name'], \
				'orien':None,\
			} ) )
			
		# add as a whole.
		self.add_nodes_from( to_add )
		


	def remove_cycles(self):
		''' removes cycles from graph '''
		
		# copy self to simple digraph.
		dig = self.di_subgraph(self.nodes())
		
		# identify cycles.
		cycles = nx.simple_cycles(dig)
		
		# loop over cycles.
		for c in cycles:
			
			# build set of edges.
			valid = dict()
			for i in range(len(c) - 1):
				if c[i] < c[i+1]:
					valid[(c[i], c[i+1])] = 0
				else:
					valid[(c[i+1], c[i])] = 0
		
			# remove least supported cycle.
			for i in range(self._edges.size):
				
				# check if part of cycle.
				if self._edges[i]['ctg_a_idx'] < self._edges[i]['ctg_b_idx']:
					key = (self._edges[i]['ctg_a_idx'], self._edges[i]['ctg_b_idx'])
				else:
					key = (self._edges[i]['ctg_b_idx'], self._edges[i]['ctg_a_idx'])
					
				if key not in valid: continue
				
				# count.
				valid[key] += 1
				
			# find least supported.
			least_a = 999999999999999
			least_b = 999999999999999
			for key in valid:
				if valid[key] < least_a:
					least_a = valid[key]
					least_b = key
			key = least_b
			
			# remove this edge.
			if self.has_edge(key[0], key[1]):
				self.remove_edge(key[0], key[1])
			else:
				self.remove_edge(key[1], key[0])
