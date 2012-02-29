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
import subprocess
from operator import itemgetter

# graph imports.
from data_structs.types import agp_dt


class ScaffAgp(nx.DiGraph):
	
	def __init__(self, agps, skip=False):
		''' Linear graph constructor.'''
		# call directed graph constructor.
		nx.DiGraph.__init__(self)
		
		# save agps.
		self.agp_list = agps
		self.agp_sz = self.agp_list.size
		
		# add nodes and edges.
		if skip == False:
			self._add_nodes()
			self._add_edges()
				

	def _add_edges(self):
		''' grabs edges from directed graph.'''
		logging.debug("adding edges.")
	
		# add edges from N components.
		for i in range(self.agp_sz):
			# skip contigs.
			if self.agp_list[i]['comp_type'] == "N":

				# sanity check.
				if i + 1 >= self.agp_list.size:
					
					# hit end, break.
					break
					
				# simplify.
				p = self.agp_list[i-1]['comp_name']
				q = self.agp_list[i+1]['comp_name']
				esize = self.agp_list[i]['comp_stop'] - self.agp_list[i]['comp_start']	
	
				# add edge.
				self.add_edge(p, q, weight=esize)
		return
			
	def _add_nodes(self):
		''' adds nodes from list.'''
		logging.debug("adding nodes.")
		
		# build list of contigs.
		nodes = []
		for i in range(self.agp_sz):
			
			# skip fragments.
			if self.agp_list[i]['comp_name'] == "fragment":
				continue
			
			# add node to graph.
			nodes.append( (self.agp_list[i]['comp_name'], {'id':i} ) )
		
		# add nodes in one shot.
		self.add_nodes_from( nodes )

