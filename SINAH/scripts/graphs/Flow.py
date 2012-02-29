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
class Flow(nx.Graph):
	
	def __init__(self, TG):
		''' creates bi-partite graph from directed one '''
		
		# call base class constructor.
		nx.Graph.__init__(self)

		# save pointer to directed graph.
		self._TG = TG

		# create arrays to organize nodes.
		self._left = list()
		self._right = list()

		# add all the nodes.
		self._add_nodes()
		
		# add the edges based on directed graph.
		self._add_edges()
		
	def pathify(self):
		''' compute maximum weight matching ''' 
		
		# call function.
		matching = nx.max_weight_matching(self, maxcardinality=True)
		
		# retrieve edges.
		elist = []
		for key in matching:
			
			# tokenize.
			if key.count("left") > 0:
				l = int(key.replace("left_",""))
				r = int(matching[key].replace("right_",""))
			else:
				r = int(key.replace("right_",""))
				l = int(matching[key].replace("left_",""))
			
			# add edge to list.
			elist.append( (r,l) )
			
		# return edgelist.
		return elist
		
	def _add_edges(self):
		''' adds edges between nodes'''	
		
		# loop over each edge.
		elist = []
		for e in self._TG.edges():
			
			# left of source.
			r = "right_%i" % e[0]
			
			# right of target.
			l = "left_%i" % e[1]
			
			# get weight.
			wt = self._TG[e[0]][e[1]]['weight']
			
			# add edge.
			elist.append( (l, r, {'weight':wt}) )
			
		# add list of edges.
		self.add_edges_from(elist)
			
		
		
	def _add_nodes(self):
		''' adds two nodes for each node in TG'''
		
		# loop over nodes.
		for n in self._TG.nodes():
			
			# add left / right.
			self._left.append( "left_%i" % n )
			self._right.append( "right_%i" % n )
			
		# add nodes.
		self.add_nodes_from(self._left)
		self.add_nodes_from(self._right)
			
