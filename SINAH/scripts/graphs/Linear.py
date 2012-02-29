'''
builds a linear graph from solution.
'''
# program imports.
from data_structs.agps import agp_dt
from graphs.MstGraph import MstGraph

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
class Linear(nx.DiGraph):
	
	def __init__(self, nodes, edges):
		
		# call base class constructor.
		nx.DiGraph.__init__(self)
		
		# save pointers to important structures.
		self._nodes = nodes
		self._edges = edges
		
		# add all the nodes right away.
		self._add_nodes()
		
	def _is_linear(self, subg=False):
		''' checks if a graph is linear '''
		
		# set active graph.
		if subg == False:
			ag = self
		else:
			ag = subg
			
		# loop over each node.
		for n in ag.nodes():
			
			# ensure its in degree < 2
			if ag.in_degree(n) >= 2:
				return False
				
			if ag.out_degree(n) >= 2:
				return False
				
		return True		
		
	def prune_negative(self):
		''' removes negative edges '''
		
		# loop over edges.
		elist = []
		for e in self.edges():
			
			# check distance.
			if self[e[0]][e[1]]['dist'] < (-1 * 3 * 250):
				elist.append(e)
				
		# remove negative.
		self.remove_edges_from(elist)
		
		
	def set_orientation(self, sol):
		''' sets orientation of nodes '''
		
		# get nodes and bundles.
		sol_nodes = sol.get_nodes()
		
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
		
	def verify(self):
		''' ensures graph is linear '''
		
		# ensure its a DAG.
		if nx.is_directed_acyclic_graph(self) == False:
			logging.error("graph is not acyclic")
			sys.exit(1)
		
		# loop over each node.
		for n in self.nodes():
			
			# ensure its in degree < 2
			if self.in_degree(n) >= 2:
				logging.error("not linear: in")
				print n, self.predecessors(n)
				sys.exit(1)
				
			if self.out_degree(n) >= 2:
				logging.error("not linear: out")
				print n, self.successors(n)
				sys.exit(1)	
				
			# ensure its orientation is set.
			if self.node[n]['orien'] != 0 and self.node[n]['orien'] != 1:
				logging.error("node orientation not set")
				sys.exit(1)
				
		# ensure the distances are set.
		for e in self.edges():
			
			# check dist.
			if 'dist' not in self[e[0]][e[1]]:
				logging.error("dist was never even set at all")
				print e
				sys.exit(1)
			if self[e[0]][e[1]]['dist'] == None:
				logging.error("distance not set: %s %s %s" % (e[0], e[1], self[e[0]][e[1]]['dist']))
				sys.exit(1)
		
	def set_gaps(self, gaps):
		''' sets gaps given dictionary '''
		
		# loop over gaps.
		for key in gaps:
			
			# assign estimate.
			if self.has_edge(key[0], key[1]) == True:
				self[key[0]][key[1]]['dist'] = gaps[key]
				
			elif self.has_edge(key[1], key[0]) == True:
				self[key[1]][key[0]]['dist'] = gaps[key]
				
			#else:
			#	logging.error("gap not in graph")
			#	sys.exit(1)
				
	def set_gaps_agp(self, gaps):
		''' sets gaps given agp '''
		
		# loop over gaps.
		for key in agp_gaps:
			
			# assign estimate.
			if self.has_edge(key[0], key[1]) == True:
				self[key[0]][key[1]]['dist'] = gaps[key]
				
			elif self.has_edge(key[1], key[0]) == True:
				self[key[1]][key[0]]['dist'] = gaps[key]
				
			else:
				logging.error("gap not in graph")
				sys.exit(1)
								
	def create_elist(self, elist):
		''' adds edges from a list '''
		
		# add edges.
		self.add_edges_from(elist)	

		
	def create_strict(self, sol):
		''' creates graph from strict ILP '''
		
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
			
			# add edge based on status.
			if sol_bundles[i]['X'] == 1:
				to_add.append( (idxa, idxb, {\
					'id':i,\
					'dist':None\
				}) )
			elif sol_bundles[i]['Y'] == 1:
				to_add.append( (idxb, idxa, {\
					'id':i,\
					'dist':None\
				}) )
		
		# add the edges at once.
		self.add_edges_from(to_add)
		
	def create_weak(self, sol):
		''' creates graph from strict ILP '''
		
		# get nodes and bundles.
		sol_nodes = sol.get_nodes()
		sol_bundles = sol.get_bundles()
		
		# add the orientation to nodes.
		for i in range(sol_nodes.size):
			
			# simplify.
			idx = sol_nodes[i]['idx']
			orien = sol_nodes[i]['orien']

			# set value.
			self.node[idx]['orien'] = orien
			
		
		# add additional computed edges.
		edges = []
		for edge in sol_bundles:
			
			# get index.
			idxa = edge['idxa']
			idxb = edge['idxb']

			# get node information.
			ctg_a_o = self.node[idxa]['orien']
			ctg_b_o = self.node[idxb]['orien']
			ctg_a_sz = self.node[idxa]['size']
			ctg_b_sz = self.node[idxb]['size']
			
			# get edge information.
			stateA = edge['A']
			stateB = edge['B']
			stateC = edge['C']
			stateD = edge['D']
			
			# check not set.
			if stateA == -1 or stateB == -1 or stateC == -1 or stateD == -1:
				logging.error("State not set for the edge\n%s", '\t'.join([str(x) for x in edge]))
				sys.exit(1)
			
			# check if there is an edge.
			edge = [-1,-1]
			if stateA == 1:
				if ctg_a_o == 0 and ctg_b_o == 0:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 1 and ctg_b_o == 1:
					edge[0] = idxb
					edge[1] = idxa
				state = 0
			elif stateD == 1:
				if ctg_a_o == 1 and ctg_b_o == 1:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 0 and ctg_b_o == 0:
					edge[0] = idxb
					edge[1] = idxa
				state = 3
			elif stateB == 1:
				if ctg_a_o == 0 and ctg_b_o == 1:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 1 and ctg_b_o == 0:
					edge[0] = idxb
					edge[1] = idxa
				state = 1
			elif stateC == 1:
				#print "debug2", ctg_a_o, ctg_b_o
				if ctg_a_o == 1 and ctg_b_o == 0:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 0 and ctg_b_o == 1:
					edge[0] = idxb
					edge[1] = idxa
				state = 2

			
			# check if invalid.
			if edge[0] == -1 or edge[1] == -1:
				continue
			
			# check if in there.
			if self.has_edge( edge[0], edge[1] ) == True:
				continue
			
			# add edge accordingly.
			edges.append( (edge[0], edge[1], {"special":True}) )
			
		# add edges from list.
		logging.debug("adding %i edges to path" % len(edges))
		self.add_edges_from(edges)
		
		# make a path.
		#self._pathify()
		
	def expand_weak(self, sol):
		''' adds special edges not part of the computed path '''
		
		# get nodes and bundles.
		sol_nodes = sol.get_nodes()
		sol_bundles = sol.get_bundles()
		
		# add additional computed edges.
		edges = []
		for edge in sol_bundles:
			
			# get index.
			idxa = edge['idxa']
			idxb = edge['idxb']

			# get node information.
			ctg_a_o = self.node[idxa]['orien']
			ctg_b_o = self.node[idxb]['orien']
			ctg_a_sz = self.node[idxa]['size']
			ctg_b_sz = self.node[idxb]['size']
			
			# get edge information.
			stateA = edge['A']
			stateB = edge['B']
			stateC = edge['C']
			stateD = edge['D']
			
			# check not set.
			if stateA == -1 or stateB == -1 or stateC == -1 or stateD == -1:
				logging.error("State not set for the edge\n%s", '\t'.join([str(x) for x in edge]))
				sys.exit(1)
			
			# check if there is an edge.
			edge = [-1,-1]
			if stateA == 1:
				if ctg_a_o == 0 and ctg_b_o == 0:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 1 and ctg_b_o == 1:
					edge[0] = idxb
					edge[1] = idxa
				state = 0
			elif stateD == 1:
				if ctg_a_o == 1 and ctg_b_o == 1:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 0 and ctg_b_o == 0:
					edge[0] = idxb
					edge[1] = idxa
				state = 3
			elif stateB == 1:
				if ctg_a_o == 0 and ctg_b_o == 1:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 1 and ctg_b_o == 0:
					edge[0] = idxb
					edge[1] = idxa
				state = 1
			elif stateC == 1:
				#print "debug2", ctg_a_o, ctg_b_o
				if ctg_a_o == 1 and ctg_b_o == 0:
					edge[0] = idxa
					edge[1] = idxb
				elif ctg_a_o == 0 and ctg_b_o == 1:
					edge[0] = idxb
					edge[1] = idxa
				state = 2

			
			# check if invalid.
			if edge[0] == -1 or edge[1] == -1:
				continue
			
			# check if in there.
			if self.has_edge( edge[0], edge[1] ) == True:
				continue
			
			# add edge accordingly.
			edges.append( (edge[0], edge[1], {"special":True}) )
			
		# add edges from list.
		logging.debug("adding %i edges to path" % len(edges))
		self.add_edges_from(edges)

		
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
		
	def order_bundles(self, bundles, sol_bundles):
		''' order all bundles by support '''
		
		# build bundle lookup.
		self._blookup = dict()
		for i in range(bundles.size):
			
			# make key.
			if bundles[i]['ctg_a_idx'] < bundles[i]['ctg_b_idx']:
				key = (bundles[i]['ctg_a_idx'], bundles[i]['ctg_b_idx'])
			else:
				key = (bundles[i]['ctg_b_idx'], bundles[i]['ctg_a_idx'])
				
			# add to lookup.
			self._blookup[key] = i
	
		# allocate weight array.
		supp = np.zeros(sol_bundles.size, dtype=np.double)
		
		# populate weight.
		keys = []
		for i in range(sol_bundles.size):
			
			# make key.
			idxa = sol_bundles[i]['idxa']
			idxb = sol_bundles[i]['idxb']
			if idxa < idxb:
				key = (idxa, idxb)
			else:
				key = (idxb, idxa)
			
			# get index.
			if key in self._blookup:
				idx = self._blookup[key]
			else:
				supp[i] = 999999
				continue
			
			# get weight.
			wt = -1
			for x in ['A', 'B', 'C', 'D']:
				if sol_bundles[i][x] == 1:
					# grab weight.
					wt = bundles[idx]["WT_%s" % x]
			if wt == -1: continue
			
			# add to supp.
			supp[i] = wt
			
			# note we used key.
			keys.append(key)
			
		# perform an indirect sort.
		abc = np.argsort(supp)
		
		# shuffle the array.
		bcopy = sol_bundles.copy()
		bcopy = bcopy[abc]
		
		# index where each key is in the argsort.
		self._bidx = dict()
		for i in range(bcopy.size):

			# make key.
			idxa = bcopy[i]['idxa']
			idxb = bcopy[i]['idxb']
			if idxa < idxb:
				key = (idxa, idxb)
			else:
				key = (idxb, idxa)
				
			# note it.
			self._bidx[key] = i
		

	def remove_cycles(self):
		''' removes cycles from graph '''	

		# identify strongly connected components.
		comps = nx.strongly_connected_components(self.di_subgraph(self.nodes()))
			
		# iterate over components.
		for comp in comps:
			
			
			# skip small.
			if len(comp) < 2: continue
			
			# get edges.
			edges = self.edges(nbunch=comp)
			
			self.remove_edges_from(edges)
			continue
			
			# make subgraph.
			g = nx.DiGraph()
			g.add_edges_from(edges)
			
			# get cycles.
			cycles = nx.simple_cycles(g)
			
			#print cycles
			
			# break each cyce.
			for cycle in cycles:
				
				# build edge list.
				elist = []
				for i in range(0, len(cycle)-1):
					if self.has_edge(cycle[i], cycle[i+1]) == True:
						elist.append( (cycle[i], cycle[i+1]) )
					else:
						elist.append( (cycle[i+1], cycle[i]) )
			
				# identify smallest.
				key = None
				val = 9999999
				for e in elist:
					if self._bidx[e] < val:
						key = e
						val = self._bidx[e]
			
				# skip.
				if key == None: continue
			
				# remove smallest.
				#print key
				self.remove_edge(key[0], key[1])
			

				
		
		
		
	def write_agp(self, agp_file):
		''' writes graph to an AGP file '''
		
		# grab some info.
		nsz = self.number_of_nodes()
		esz = self.number_of_edges()
		tsz = nsz + esz
		
		# create AGP array.
		self.agp = np.zeros(tsz, dtype=agp_dt)
		
		# locate connected components.
		components = nx.weakly_connected_components(self)
		
		# loop over components.
		idx = 0
		agp_idx = 0
		for component in components:
			
			# create subgraph.
			subg = nx.DiGraph()
			subg.add_nodes_from(component)
			
			# add edges.
			for edge in self.edges(component):
				subg.add_edge(edge[0],edge[1],dist=self[edge[0]][edge[1]]['dist'])
			
			# compute topological order.
			path = nx.topological_sort(subg)
			
			# validate path.
			pn = path[0]
			for i in range(1,len(path)):
				if subg.has_edge(pn, path[i]) == False:
					logging.error("Not a linear graph.")
					sys.exit()
				pn = path[i]
				
			
			# setup scaffold vars.
			scaf_idx = 1
			scaf_start = 1
			scaf_stop = self.node[path[0]]['size']
			
			# add first node.
			self.agp[agp_idx]['scaf_name'] = "scaf_%i" % idx
			self.agp[agp_idx]['scaf_start'] = scaf_start
			self.agp[agp_idx]['scaf_stop'] = scaf_stop
			self.agp[agp_idx]['scaf_idx'] = scaf_idx
			self.agp[agp_idx]['comp_type'] = "W"
			self.agp[agp_idx]['comp_name'] = self.node[path[0]]['name']
			self.agp[agp_idx]['comp_start'] = 1
			self.agp[agp_idx]['comp_stop'] = self.node[path[0]]['size']
			self.agp[agp_idx]['comp_orien'] = self.node[path[0]]['orien']
			self.agp[agp_idx]['comp_linkage'] = 0
			
			# increment pointers.
			agp_idx += 1
			pn = path[0]
			scaf_idx += 1
			scaf_start = scaf_stop + 1
			
			# fill AGP.
			for i in range(1,len(path)):

				# defualt gapsize to 7.
				gapsz = subg[pn][path[i]]['dist']
				if gapsz < 0:
					gapsz = 7

				# get stop location.
				scaf_stop = scaf_start + gapsz
				
				# add gap info.
				self.agp[agp_idx]['scaf_name'] = "scaf_%i" % idx
				self.agp[agp_idx]['scaf_start'] = scaf_start
				self.agp[agp_idx]['scaf_stop'] = scaf_stop
				self.agp[agp_idx]['scaf_idx'] = scaf_idx
				self.agp[agp_idx]['comp_type'] = "N"
				self.agp[agp_idx]['comp_name'] = "fragment"
				self.agp[agp_idx]['comp_start'] = 1
				self.agp[agp_idx]['comp_stop'] = gapsz
				self.agp[agp_idx]['comp_linkage'] = 0
				
				# increment pointers.
				agp_idx += 1
				scaf_idx += 1
				scaf_start = scaf_stop + 1
				
				# set stop.
				scaf_stop = scaf_start + self.node[path[i]]['size']
				
				# add nodes.
				self.agp[agp_idx]['scaf_name'] = "scaf_%i" % idx
				self.agp[agp_idx]['scaf_start'] = scaf_start
				self.agp[agp_idx]['scaf_stop'] = scaf_stop
				self.agp[agp_idx]['scaf_idx'] = scaf_idx
				self.agp[agp_idx]['comp_type'] = "W"
				self.agp[agp_idx]['comp_name'] = self.node[path[i]]['name']
				self.agp[agp_idx]['comp_start'] = 1
				self.agp[agp_idx]['comp_stop'] = self.node[path[i]]['size']
				self.agp[agp_idx]['comp_orien'] = self.node[path[i]]['orien']
				
				# increment pointers.
				pn = path[i]
				agp_idx += 1
				scaf_idx += 1
				scaf_start = scaf_stop + 1

		
			# increment scaffold idx
			idx += 1
			
		# write to file.
		fout = open(agp_file, "w")
		
		# write each entry.
		z = len(agp_dt.names)
		for i in range(self.agp.size):
			
			# format result.
			tmp = self.agp[i]
			if tmp['comp_type'] == "W":
				# get orientation.
				if tmp["comp_orien"] == 0:
					o = "+"
				else:
					o = "-"
					
				# write contig.
				txt = str(tmp['scaf_name']) + "\t"
				txt += str(tmp['scaf_start']) + "\t"
				txt += str(tmp['scaf_stop']) + "\t"
				txt += str(tmp['scaf_idx']) + "\t"
				txt += str(tmp['comp_type']) + "\t"
				txt += str(tmp['comp_name']) + "\t"
				txt += str(tmp['comp_start']) + "\t"
				txt += str(tmp['comp_stop']) + "\t"
				txt += o + "\n"
				
			else:
				# get linkage.
				if tmp['comp_linkage'] == 0:
					o = "no"
				else:
					o = "yes"
				
				# write gap.
				txt = str(tmp['scaf_name']) + "\t"
				txt += str(tmp['scaf_start']) + "\t"
				txt += str(tmp['scaf_stop']) + "\t"
				txt += str(tmp['scaf_idx']) + "\t"
				txt += str(tmp['comp_type']) + "\t"
				txt += str(tmp['comp_stop'] - tmp['comp_start']) + "\t"
				txt += str(tmp['comp_name']) + "\t"
				txt += o + "\n"
								
			fout.write(txt)
			
		# close file.
		fout.close()

	def di_subgraph(self, nbunch):
		"""Return the subgraph induced on nodes in nbunch."""
		
		# create new subgraph.
		g = nx.DiGraph()
		
		# add edges and nodes.
		g.add_edges_from( self.edges(nbunch=nbunch) )
		
		return g

	def pathify(self):
		''' makes a path out of graph '''
		
		
		# build stack of unverified components.
		stack = nx.weakly_connected_components(self)
		
		# loop until stack is empty.
		while len(stack) > 0:
			
			print len(stack)
			
			# pop off stack.
			comp = stack.pop()

			# skip singeltons.
			if len(comp) < 3:
				continue
				
			# create proper subgraph.
			subg1 = self.di_subgraph(nbunch=comp)
		
			# make input for mst graph.
			successors = nx.to_dict_of_lists(subg1) 
			
			# make the MstGraph.
			mst_graph = MstGraph(successors)
			
			# calculate the MST.
			mst = mst_graph.mst()
			
			# build graph from MST.
			subg2 = nx.DiGraph()
			for edge in mst.iteredges():
				subg2.add_edge(edge[0], edge[1], {'weight':-1})
				
			# ensure its acylclic.
			if nx.is_directed_acyclic_graph(subg2) == False:
				logging.error("pathify: cycle in graph")
				sys.exit(1)
				
			# use dikstra all pairs to find backbone.
			paths = nx.all_pairs_dijkstra_path(subg2, weight='weight')
			
			# determine longest.
			longest_a = -1
			longest_b = -1
			for path in paths:
				if len(paths[path]) > longest_a:
					longest_a = len(paths[path])
					longest_b = path
			path = paths[longest_b]
			
			# remove edges 1 adjacent to path.
			to_remove = []
			for e in subg2.edges():
				test1 = e[0] in path
				test2 = e[1] in path
				if test1 + test2 == 1:
					to_remove.append(e)
					print "removing"
					
			if self._is_linear(subg=subg2) == False:
				
				print 
				sys.exit()
					
			# remove edges.
			subg2.remove_edges_from(to_remove)
			
			# calculate components.
			comps = nx.weakly_connected_components(subg1)
			
			# add unverified to stack.
			for comp in comps:
				
				# build new graph.
				subg2 = self.di_subgraph(nbunch=comp)
				
				# verify.
				if self._is_linear(subg=subg2) == False:
					stack.append(comp)
	

def bellman_ford(G, source, weight = 'weight'):
	if source not in G:
		raise KeyError("Node %s is not found in the graph"%source)
	numb_nodes = len(G)
	dist = {source: 0}
	pred = {source: None}
	if numb_nodes == 1:
	   return pred, dist
   
	if G.is_multigraph():
		def get_weight(edge_dict):
			return min([eattr.get(weight,1) for eattr in edge_dict.values()])
	else:
		def get_weight(edge_dict):
			return edge_dict.get(weight,1)
	for i in range(numb_nodes):
		no_changes=True
		# Only need edges from nodes in dist b/c all others have dist==inf
		for u, dist_u in list(dist.items()): # get all edges from nodes in dist
			for v, edict in G[u].items():  # double loop handles undirected too
				dist_v = dist_u + get_weight(edict)
				if v not in dist or dist[v] > dist_v:
					dist[v] = dist_v
					pred[v] = u
					no_changes = False
		if no_changes:
			break
	else:
		raise nx.NetworkXUnbounded("Negative cost cycle detected.")
	return pred, dist
