#!/usr/bin/python
'''
Builds a graph from node and edge file.
'''
# system imports.
import sys
import logging
import copy
import numpy as np
import networkx as nx
import matplotlib.pylab as plt
import itertools
import math

class ScaffBundle(nx.Graph):
	
	def __init__(self):
		''' extends normal graph class with more methods'''
		
		# call multi graph constructor.
		nx.Graph.__init__(self)
		
	
	def bootstrap(self, node_list, edge_list):
		''' adds nodes and edges to the bundle graph '''
	
		# add nodes and edges to self.
		self._add_nodes(node_list)
		self._add_edges(edge_list)
		
	def induced_bigraph(self, nbunch, x):
		''' creates induced biconnected graph, x is the cut node
		and only edges between x and the subgraph are added.
		 '''
		 
		# ensure x is not in nbunch.
		if x in nbunch:
			nbunch.remove(x)
			 
		# create subgraph.
		subg = self.subgraph(nbunch)
		
		# add x.
		subg.add_node(x)
		
		# add edges from x into subgraph.
		for y in self.neighbors(x):
			
			# check if neighbor in set.
			if y in nbunch:
				subg.add_edge(x, y, {'id':self[x][y]['id']})
				
		# return subgraph.
		return subg
		
	def induced_trigraph(self, nbunch, x, y):
		''' creates induced triconnected graph, x,y are the cut node
		and only edges between x or y and the subgraph are added.
		 '''
		 
		# ensure x is not in nbunch.
		if x in nbunch:
			nbunch.remove(x)
		if y in nbunch:
			nbunch.remove(y)
			 
		# create subgraph.
		subg = self.subgraph(nbunch)
		
		# add x, y.
		subg.add_node(x)
		subg.add_node(y)
		
		# add edges from x into subgraph.
		for z in self.neighbors(x):
			
			# check if neighbor in set.
			if z in nbunch:
				subg.add_edge(x, z, {'id':self[x][z]['id']})
				
		# add edges from y into subgraph.
		for z in self.neighbors(y):
			
			# check if neighbor in set.
			if z in nbunch:
				subg.add_edge(y, z, {'id':self[y][z]['id']})
				
		# validate subgraph.
		for x in subg.edges():
			subg[x[0]][x[1]]['id']
				
		# return subgraph.
		return subg


	def get_triangles(self):
		''' gets cycles in graph as a whole. '''
		
		# neighbor set.
		neibs = {}
		
		# set node list
		node_list = self.nodes()
		
		# loop over each node.
		for p in node_list:
			# get neighbors.
			tmp = set()
			for q in self[p]:
				
				# add to set.
				tmp.add(q)
				
			# add set to dict.
			neibs[p] = tmp
		
		# find all triangles.
		tris = []
		for p in node_list:
			
			# look for intersection in any two neighbors.
			for pair in itertools.permutations(neibs[p],2):
				
				# simplify.
				x = pair[0]
				y = pair[1]
				
				setx = neibs[pair[0]]
				sety = neibs[pair[1]]
				
				# check if neighbor in set.
				if y in setx:
					# found cycle.
					tris.append([p, x, y, p])
					
		# return cycles.
		return tris
			

	def _add_nodes(self, node_list):
		''' adds nodes from list.'''
					
		# build nodes
		to_add = list(node_list[:]['node_idx'])
		
		# add in set.
		self.add_nodes_from(to_add)

		
	def _add_edges(self, bundle_list):
		''' adds edges from list.'''
		
		# build list to add.
		toadd = []
		
		# loop over each bundle.
		for i in range(bundle_list.size):
			
			# add edge and to list.
			toadd.append( (bundle_list[i]['ctg_a_idx'], bundle_list[i]['ctg_b_idx'], {'id':i, 'a_idx':bundle_list[i]['ctg_a_idx']}) )
			
			
		# add merged edges into graph.
		self.add_edges_from(toadd)

	def draw(self, name, subg=False):
		''' draw an interactive graph.'''
		# write dot file.
		dot = "./tmp.dot"
		pic = "/home/jlindsay/central/transfer/%s.jpg" % name
		
		if subg == False:
			nx.write_dot(self, dot)
		else:
			nx.write_dot(subg, dot)
		
		# convert to picture using neato.
		subprocess.call(["neato", "-Tjpg", "-o", pic, dot])
		
		# remove dot.
		subprocess.call(["rm","-f",dot])
