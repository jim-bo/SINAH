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

# program imports.
from graphs.MstGraph import MstGraph
from data_structs.types import agp_dt
from data_structs.types import edge_dt


############ functions #########

def determine_dist(aid, bid, state, edge, node1, node2):
	''' given a state, determines distance between reads implied by that state	'''
	
	# determine column to work with.
	if aid < bid:
		
		# calc distance based on state.
		if state == 0:
			d1 = node1['ctg_width'] - edge['read_a_right_pos']
			d2 = edge['read_b_left_pos']
		elif state == 1:
			d1 = node1['ctg_width'] - edge['read_a_right_pos']
			d2 = node2['ctg_width'] - edge['read_b_right_pos']
		elif state == 2:
			d1 = edge['read_a_left_pos']
			d2 = edge['read_b_left_pos']
		else:
			d1 = edge['read_a_left_pos']
			d2 = node2['ctg_width'] - edge['read_b_right_pos']
			
	else:
		
		# calc distance based on state.
		if state == 0:
			d1 = edge['read_a_left_pos']
			d2 = node2['ctg_width'] - edge['read_b_right_pos']
		elif state == 1:
			d1 = node1['ctg_width'] - edge['read_a_right_pos']
			d2 = node2['ctg_width'] - edge['read_b_right_pos']
		elif state == 2:
			d1 = edge['read_a_left_pos']
			d2 = edge['read_b_left_pos']
		else:
			d1 = node1['ctg_width'] - edge['read_a_right_pos']
			d2 = edge['read_b_left_pos']
			
	# sanity check.
	if d1 < 0:
		print "bad dist: d1"
		print d1
		print state
		print "bad dist: d1"
		sys.exit(1)
	if d2 < 0:
		print "bad dist: d2"
		print d2
		print state
		print aid, bid
		print "bad dist: d2"
		sys.exit(1)
		
	# finalize dist.
	dist = edge['insert_size'] - d1 - d2
		
	# return dist.
	return d1, d2
	
def determine_state(aid, bid, ao, bo):
	''' determines class of edges '''
	
	# check if type 1
	if aid < bid:
			
		# determine happyness class.
		if ao == 0 and bo == 0:
			type = 0
		elif ao == 0 and bo == 1:
			type = 1
		elif ao == 1 and bo == 0:
			type = 2
		elif ao == 1 and bo == 1:
			type = 3
			
	# no its type two.
	else:
		
		# determine happyness class.
		if ao == 1 and bo == 1:
			type = 0
		elif ao == 0 and bo == 1:
			type = 1
		elif ao == 1 and bo == 0:
			type = 2
		elif ao == 0 and bo == 0:
			type = 3
	
	# return type.
	return type

def update_pos(node, edge):
	''' updates the edge Lpos, Rpos based on flip status of node '''
	
	# create new edge.
	tmpe = np.zeros(1, dtype=edge_dt)[0]
	
	# copy data into it.
	for x in edge_dt.fiels:
		tmpe[x] = edge[e]
	
	# check orientation status.
	if node['ctg_orien'] == 0:
		return tmpe
		
	# flip them.
	

############ class ###########
class ScaffLinear(nx.DiGraph):
	
	def __init__(self, nodes, bundles, edges, nobreak=False, pseudo=False):
		''' Linear graph constructor.'''
		# call directed graph constructor.
		nx.DiGraph.__init__(self)
		
		# skip if its a pseudo linear graph (usually made from output of other program)
		if pseudo == True:
			return
		
		# add pointers.
		self.node_list = nodes
		self.bundle_list = bundles
		
		# node idx.
		self.new_idx = 0
		
		# add bundle list.
		self.node_sz = self.node_list.size
		self.bundle_sz = self.bundle_list.size
		
		# add edge list.
		self.edge_list = edges
		self.edge_sz = self.edge_list.size
		
		# picture tracker.
		self.pic_idx = 0
			
		# compute gap estimates.
		self.gaps = self._gap_avg(2000)
		#self.gaps = self._gap_EM()
				
		# add nodes and edges.
		self._add_nodes()
		self._add_edges(self.gaps)
			
		# choose a path.
		self._pathify()
			
		# verify its linearity.
		for comp in nx.weakly_connected_components(self):
			
			# make subgraph.
			subg = self.di_subgraph(comp)
			
			# check linear,.
			if self._is_linear(subg) == False:
				self.draw("failed",subg=subg)
				logging.error("non linear path")
				sys.exit()
		
		# check all nodes are present.
		for n in range(self.node_sz):
			if self.has_node(n) == False:
				logging.error("Missing node in linear graph.\n%s\n" % str(self.node_list[n]))
				sys.exit(1)


	def _greedy(subg):
		''' return digraph of max scoring edge for each node'''
		# create new graph.
		g = nx.DiGraph()
		
		# loop over each each node.
		for n in subg.nodes():
			
			# find max score incoming link.
			max_score = -99999999999
			max_p = None
			
			# loop over predecessors.
			for p in subg.predecessors(n):
				
				# score.
				if subg[p][n]['weight'] > max_score:
					max_score = subg[p][n]['weight']
					max_p = p
					
			# add the max edge.
			g.add_edge(max_p, n, {'weight':max_score})
			
		# return graph.
		return g
			
			
	def _contract(self, subg, cycle):
		'''Given a cycle in our graph, contract it into a single node.

		Returns a tuple (id, graph). The graph is a new Digraph instance
		containing no nodes from the cycle, with one extra new node created to
		represent the cycle. The id is the id of the new node.
		'''
		
		# create the new graph.
		g = nx.DiGraph()
		
		# add all nodes not in cycle.
		for n in subg.nodes():
			g.add_node(n)
			
		# create cycle set.
		cset = set(cycle)
			
		# add all edges not incident to cycle.
		for edge in subg.edges():
			
			# skip incidents.
			if edge[0] in cset or edge[1] in cset:
				continue
				
			# add edge.
			g.add_edge(edge)
			
		# determine max pred weight of cycle.
		max_pred = - 999999999999
		for n in cycle:
			for p in subg.predecessors(n):
				if subg[p][n]['weight'] > max_score:
					max_pred = subg[p][n]['weight']
					
		# determin max succ weight.
		max_succ = - 999999999999
		for n in cycle:
			for p in subg.successors(n):
				if subg[p][n]['weight'] > max_score:
					max_succ = subg[p][n]['weight']
					
		# create new node.
		x = "fake_%i" % self.new_idx
		self.new_idx += 1
		
		# add edges to new node.
		for n in cycle:
			
			# add incomming edges.
			for p in subg.predecessors(n):
				g.add_edge(p, x, {'weight':max_pred})
			
			# add outgoing edges.
			for p in subg.successors(n):
				g.add_edge(x, p, {'weight':max_succ})
		
		# return the graph.
		return x, subg.edges(), g
		
	def dmst(self, g):
		'''Return the MST of this Digraph using the Chu-Liu-Edmonds algorithm'''
		
		# find max scoring digraph.
		candidate = self._greedy(g)
		
		# find cycle in this graph.
		cycles = nx.simple_cycles(candidate)
		
		# if no cycle return graph.
		if cycles == []:
			return candidate
		
		# contract this graph.
		new_id, old_edges, compact = self._contract(cycle)
		
		# merge the contracted graph recursively.
		merged = self.merge(self.dmst(compact), new_id, old_edges, cycle)
		return merged



	def _pathify(self):
		''' makes a path out of graph '''
		
		
		# identify paths.
		for comp in nx.weakly_connected_components(self):
			
			# skip singeltons.
			if len(comp) == 1:
				continue
		
			# keep track to edges to remove.
			to_remove = []
			to_add = []
			
			# create proper subgraph.
			subg1 = self.di_subgraph(nbunch=comp)
		
			# make input for mst graph.
			successors = nx.to_dict_of_lists(subg1) 
			
			# make the MstGraph.
			mst_graph = MstGraph(successors)
			
			# calculate the MST.
			mst = mst_graph.mst()
			
			# build our graph.
			subg2 = nx.DiGraph()
			
			for edge in mst.iteredges():
				subg2.add_edge(edge[0], edge[1])
			
			# get subgraph edges.
			e_edges = set(subg1.edges())		
			n_edges = set(subg2.edges())
			
			# remove edges not in subg2.
			to_remove += e_edges.difference(n_edges)
						
			# check if reduced graph is non-linear.
			if self._is_linear(subg2) == False:
						
				# debug stuff.
				x = subg2.number_of_nodes()
				if x > 20:
					logging.debug("linearizing a component of size %i" % x)
						
				# compute comps in MST.
				for comp3 in nx.weakly_connected_components(subg2):
					
					# make subgraph.
					subg3 = subg2.subgraph(nbunch=comp3)
					
					# sanity check its a DAG.
					if nx.is_directed_acyclic_graph(subg3) == False:
						logging.error("MST component has a cycle?")
						sys.exit(1)

									
					# compute the topological order
					root = nx.topological_sort(subg3)[0]
					
					# perform a DFS from root to see if we can merge paths.
					stack = [root]
					visited = {}
					
					# loop until stack is empty.
					while stack != []:
						
						# peek at node from stack.
						n = stack[-1]
						
						# push all its successors onto stack.
						visits = 0
						for p in subg3.successors(n):
							
							# check if already visited.
							if p not in visited:
								stack.append(p)
								visits += 1
								
						# if all successors visited process.
						if visits == 0:
							
							# pop node from stack.
							n = stack.pop()
							
							# get its successors.
							succs = subg3.successors(n)
							
							# store contig width if last.
							if len(succs) == 0:
								visited[n] = (float(self.node_list[n]['ctg_width']), n)
								continue
								
							# store width of successor, gap size and ctg width if only 1 succ.
							if len(succs) == 1:
								tmp = float(visited[p][0]) + float(self.gaps[n][p]) + float(self.node_list[n]['ctg_width'])
								visited[n] = (tmp, visited[p][1])
								continue
							
							
							# find the longest initial gap.
							gap_sizes = sorted([(p, self.gaps[n][p]) for p in succs], key=itemgetter)
							max_p = gap_sizes[-1][0]
							max_gap = gap_sizes[-1][1]
							
							# check if we can fit all successors within this gap.
							pos = 0.0
							for entry in gap_sizes:
								
								# simplify.
								p = entry[0]
								gap_sz = entry[1]
								
								# skip self.
								if p == max_p:
									continue
									
								# increment count.
								pos += float(visited[p][0]) + float(gap_sz)
								
							# always remove all edges.
							for p in succs:
								to_remove.append( (n,p) )
							
							# check if final position is less than max_gap.
							if pos < (max_gap + ( 3 * 250)):
								
								# merge the components.
								pnode = n
								p_pos = 0
								for entry in gap_sizes:
								
									# add edge.
									gap_sz = self.gaps[n][entry[0]] - p_pos
									to_add.append( (pnode, entry[0], {'gap_size': gap_sz}) )
									
									# add gap size.
									if pnode not in self.gaps and entry[0] not in self.gaps[pnode]:
										self.gaps[pnode][entry[0]] = gap_sz
										self.gaps[entry[0]][pnode] = gap_sz
									
									# update pnode.
									pnode = visited[entry[0]][1]
									
									# remove edge.
									to_remove.append( (n, entry[0]) )
										
								# set size to be the largest gap size plus rest of that.
								tmp = float(self.node_list[n]['ctg_width']) + max_gap + visited[max_p][1]
								visited[n] = (tmp, n)
														
							else:
															
								# set size to be ctg width.
								visited[n] = (float(self.node_list[n]['ctg_width']), n)


			# test subg2.
			self.draw("subg2_pre",subg2)
			subg2.remove_edges_from(to_remove)
			subg2.add_edges_from(to_add)
			if self._is_linear(subg2) == False:
				logging.error("mistake in pathify calculation")
				self.draw("subg1",subg1)
				self.draw("subg2",subg2)
				
				print "to_remove"
				print to_remove
				
				print "to_add"
				print to_add
				
				logging.error("bad stuff: linear.py")
				sys.exit(1)

			# remove all edges.
			#logging.debug("removing bad edges")
			self.remove_edges_from(to_remove)	
			
			# add all edges.
			#logging.debug("adding computed edges")
			self.add_edges_from(to_add)

	def _gap_avg(self, inst):
		''' computes gap estimate '''
		
		# store gaps in this dictionary.
		gaps = {}
		
		# identifies edges for each bundle.
		bunchs = {}
		
		# loop over edges.
		for i in range(self.edge_sz):
			
			# skip invalid.
			if self.edge_list[i]['invalid'] == True:
				continue
			
			# make key.
			key = tuple(sorted([self.edge_list[i]['ctg_a_idx'], self.edge_list[i]['ctg_b_idx']]))
			
			# add count.
			if key not in bunchs:
				bunchs[key] = []
			bunchs[key].append(i)
			
		# loop over edge pair.
		for bundle in self.bundle_list:
			
			# make key.
			key = tuple( sorted( [bundle['ctg_a_idx'], bundle['ctg_b_idx']] ) )
			aid = bundle['ctg_a_idx']
			bid = bundle['ctg_b_idx']
			ctg_a_sz = self.node_list[aid]['ctg_width']
			ctg_b_sz = self.node_list[bid]['ctg_width']
			
			if key not in bunchs:
				gaps[key[0]][key[1]] = 500
				gaps[key[1]][key[0]] = 500
				continue				
			
			# loop over each edge.
			tot = 0.0
			cnt = 0.0
			for i in bunchs[key]:
				
				# simplify vars.
				edge = self.edge_list[i]
				ao = edge['read_a_orien']
				bo = edge['read_b_orien']
				
				al = edge['read_a_left_pos']
				ar = edge['read_a_right_pos']
				bl = edge['read_b_left_pos']
				br = edge['read_b_right_pos']

				# arrange contig size.
				if aid < bid:
					asz = ctg_a_sz
					bsz = ctg_b_sz
				else:
					asz = ctg_b_sz
					bsz = ctg_a_sz					

				# is this type one.
				da = db = None
				if aid < bid:
					
					# determine happyness class.
					if ao == 0 and bo == 0:
						da = asz - ar
						db = bl
					elif ao == 0 and bo == 1:
						da = asz
						db = bsz - br
					elif ao == 1 and bo == 0:
						da = al
						db = bl
					elif ao == 1 and bo == 1:
						da = al
						db = bsz - br
					
				# no its type two.
				else:
					# determine happyness class.
					if ao == 1 and bo == 1:
						da = al
						db = bsz - br
					elif ao == 0 and bo == 1:
						da = asz - ar
						db = bsz - br
					elif ao == 1 and bo == 0:
						da = al
						db = bl
					elif ao == 0 and bo == 0:
						da = asz - ar
						db = bl

				# track it.
				tot += float(inst - da - db)
				cnt += 1.0		
		
			# average it.
			avg = tot / cnt
			
			# store it.
			if key[0] not in gaps:
				gaps[key[0]] = {}
			if key[1] not in gaps:
				gaps[key[1]] = {}
				
			if key[1] not in gaps[key[0]]:
				gaps[key[0]][key[1]] = None
			if key[0] not in gaps[key[1]]:
				gaps[key[1]][key[0]] = None
				
			gaps[key[0]][key[1]] = avg
			gaps[key[1]][key[0]] = avg
			
			#print avg
				
		# return the estimates.
		return gaps	

	def _gap_EM(self):
		''' computes gap estimate '''
		
		# store gaps in this dictionary.
		gaps = {}
		
		# identifies edges for each bundle.
		bunchs = {}
		
		# loop over edges.
		for i in range(self.edge_sz):
			
			# skip invalid.
			if self.edge_list[i]['invalid'] == True:
				continue
			
			# make key.
			key = tuple(sorted([self.edge_list[i]['ctg_a_idx'], self.edge_list[i]['ctg_b_idx']]))
			
			# add count.
			if key not in bunchs:
				bunchs[key] = []
			bunchs[key].append(i)
			
		# make temp edge.
		temp_edge = np.zeros(1, dtype=edge_dt)[0]
			
		# loop over edge pair.
		for bundle in self.bundle_list:
			
			# check we have an edge.
			if bundle['X'] + bundle['Y'] < 1:
				continue
			
			# make key.
			key = tuple( sorted( [bundle['ctg_a_idx'], bundle['ctg_b_idx']] ) )
	
			# initialize u.
			u = self.edge_list[bunchs[key][0]]['insert_size']
		
			# loop over each edge.
			tot = 0.0
			cnt = 0.0
			for i in bunchs[key]:
				
				# simplify.
				edge = self.edge_list[i]
				aid = edge['ctg_a_idx']
				bid = edge['ctg_b_idx']
				node1 = self.node_list[aid]
				node2 = self.node_list[bid]
				ao = node1['ctg_orien']
				bo = node2['ctg_orien']
				
				# make updated edge.
				tmp_edge = edge
				if node1['ctg_orien'] == 1:
					tmp_edge['read_a_orien'] = 1 - tmp_edge['read_a_orien']
				if node2['ctg_orien'] == 1:
					tmp_edge['read_b_orien'] = 1 - tmp_edge['read_b_orien']
					
				
				# get state.
				state = determine_state(aid, bid, ao, bo)
				
				# get edge distance.
				d1, d2 = determine_dist(aid, bid, state, tmp_edge, node1, node2)

				# track it.
				tot += float(u - d1 - d2)
				cnt += 1.0		
		
			# average it.
			x = tot / cnt
			
			
			"shant be here"
			sys.exit()

			#print avg
				
		# return the estimates.
		return gaps	

	def di_subgraph(self, nbunch):
		"""Return the subgraph induced on nodes in nbunch."""
		
		# create new subgraph.
		g = nx.DiGraph()
		
		# add edges and nodes.
		g.add_edges_from( self.edges(nbunch=nbunch) )
		
		return g
		
		
		

	def _break_cycles(self, subg):
		''' removes back edges in cycles according to worst gap size.'''
		
		# identify cycles.
		cycles = nx.simple_cycles(subg)
		
		# loop till no more cycles.
		while cycles != []:
		
			# break all cycles somehow.
			for cycle in cycles:
				
				# warn about 3 cycles.
				if len(cycle) - 1 == 3:
					logging.warning("Found a 3 cycle.")
					
				# note.
				logging.debug("found a %i cycle, breaking it" % len(cycle))
					
				# break last edge.
				cycle.reverse()
				if subg.has_edge(cycle[0], cycle[1]):
					self.remove_edge(cycle[0], cycle[1])
					subg.remove_edge(cycle[0], cycle[1])
					
				elif subg.has_edge(cycle[1], cycle[0]):
					self.remove_edge(cycle[1], cycle[0])
					subg.remove_edge(cycle[1], cycle[0])
					
			# identify cycles.
			cycles = nx.simple_cycles(subg)
				
		# done.
		return subg

	def _add_edges(self, gaps):
		''' grabs path from vars.'''
		logging.debug("adding edges.")
		
		# loop over bundle_list.
		edges = []
		for bundle in self.bundle_list:
			
			# get index.
			ida = bundle['ctg_a_idx']
			idb = bundle['ctg_b_idx']
			
			# post filter gaps with are insane (i.e. > 3 std dev)
			if gaps[ida][idb] < -(1396 * (3 *250)) or gaps[ida][idb] > (1396 * (3 *250)):
				continue

			# get bundle.
			if bundle['X'] == 1:
				edges.append((ida, idb, {"gap_size":gaps[ida][idb]}))
			
			elif bundle['Y'] == 1:
				edges.append((idb, ida, {"gap_size":gaps[ida][idb]}))

		
		# add edges.
		self.add_edges_from(edges)

	def _add_edges_old(self):
		''' grabs path from vars.'''
		logging.debug("adding edges.")
		
		# loop over bundle_list.
		edges = []
		for edge in self.bundle_list:
			
			# get index.
			id1 = edge['ctg_a_idx']
			id2 = edge['ctg_b_idx']
			
			# calculate gapsize.
			gz = 1000

			# get node information.
			ctg_a_o = self.node_list[id1]['ctg_orien']
			ctg_b_o = self.node_list[id2]['ctg_orien']
			ctg_a_sz = self.node_list[id1]['ctg_width']
			ctg_b_sz = self.node_list[id2]['ctg_width']
			
			# get edge information.
			stateA = edge['ST_A']
			stateB = edge['ST_B']
			stateC = edge['ST_C']
			stateD = edge['ST_D']
			
			# check not set.
			if stateA == -1 or stateB == -1 or stateC == -1 or stateD == -1:
				logging.error("State not set for the edge\n%s", '\t'.join([str(x) for x in edge]))
				sys.exit(1)
			
			# check if there is an edge.
			edge = [-1,-1]
			if stateA == 1:
				if ctg_a_o == 0 and ctg_b_o == 0:
					edge[0] = id1
					edge[1] = id2
				elif ctg_a_o == 1 and ctg_b_o == 1:
					edge[0] = id2
					edge[1] = id1
				state = 0
			elif stateD == 1:
				if ctg_a_o == 1 and ctg_b_o == 1:
					edge[0] = id1
					edge[1] = id2
				elif ctg_a_o == 0 and ctg_b_o == 0:
					edge[0] = id2
					edge[1] = id1
				state = 3
			elif stateB == 1:
				if ctg_a_o == 0 and ctg_b_o == 1:
					edge[0] = id1
					edge[1] = id2
				elif ctg_a_o == 1 and ctg_b_o == 0:
					edge[0] = id2
					edge[1] = id1
				state = 1
			elif stateC == 1:
				#print "debug2", ctg_a_o, ctg_b_o
				if ctg_a_o == 1 and ctg_b_o == 0:
					edge[0] = id1
					edge[1] = id2
				elif ctg_a_o == 0 and ctg_b_o == 1:
					edge[0] = id2
					edge[1] = id1
				state = 2

			
			# check if invalid.
			if edge[0] == -1 or edge[1] == -1:
				continue
			
			# add edge accordingly.
			edges.append( (edge[0], edge[1], {"gap_size":gz}) )

		
		# add edges.
		self.add_edges_from(edges)
		

	def _is_linear(self, subg):
		''' checks if provided graph is linear'''
		
		# loop over all nodes.
		fail = True
		for n in subg.nodes():
			
			# check in degree/outdegree.
			if subg.in_degree(n) > 1:
				fail = False

			if subg.out_degree(n) > 1:
				fail = False	
				
		# return result.
		return fail
	
	def _break_nonlinear(self, subg=False):
		''' checks edges are linear.'''
		
		# ensure they are linear.
		to_break = []
		fail = False
		for n in range(self.node_sz):
			
			# check in degree/outdegree.
			if self.in_degree(n) > 1:
				logging.warning("node %i has bad in_degree" % n)
				to_break.append(n)
				fail = True

			if self.out_degree(n) > 1:
				logging.warning("node %i has bad out_degree" % n)
				to_break.append(n)
				fail = True
		
		# break connected to neighbors.
		to_remove = set()
		for n in to_break:
			for x in self.out_edges(n):
				to_remove.add(x)
			for x in self.in_edges(n):
				to_remove.add(x)
				
		# remove edges.
		self.remove_edges_from(list(to_remove))
		
		return False
				
	def _calc_gap(self, aid, bid, inst):
		''' calculates weight of edge from multigraph '''
		# get node information.
		ctg_a_o = self.node_list[aid]['ctg_orien']
		ctg_b_o = self.node_list[bid]['ctg_orien']
		ctg_a_sz = self.node_list[aid]['ctg_width']
		ctg_b_sz = self.node_list[bid]['ctg_width']
		
		# loop over each edge in multigraph.
		tot = 0.0
		cnt = 0.0
		for cid in self.scaffmulti[aid][bid]:
			# get read pair.
			idx = self.scaffmulti[aid][bid][cid]['id']
			edge = self.edge_list[idx]
			
			# simplify vars.
			ao = edge['read_a_orien']
			bo = edge['read_b_orien']
			
			al = edge['read_a_left_pos']
			ar = edge['read_a_right_pos']
			bl = edge['read_b_left_pos']
			br = edge['read_b_right_pos']

			# arrange contig size.
			if aid < bid:
				asz = ctg_a_sz
				bsz = ctg_b_sz
			else:
				asz = ctg_b_sz
				bsz = ctg_a_sz					

			# is this type one.
			da = db = None
			if aid < bid:
				
				# determine happyness class.
				if ao == 0 and bo == 0:
					da = asz - ar
					db = bl
				elif ao == 0 and bo == 1:
					da = asz
					db = bsz - br
				elif ao == 1 and bo == 0:
					da = al
					db = bl
				elif ao == 1 and bo == 1:
					da = al
					db = bsz - br
				
			# no its type two.
			else:
				# determine happyness class.
				if ao == 1 and bo == 1:
					da = al
					db = bsz - br
				elif ao == 0 and bo == 1:
					da = asz - ar
					db = bsz - br
				elif ao == 1 and bo == 0:
					da = al
					db = bl
				elif ao == 0 and bo == 0:
					da = asz - ar
					db = bl

			# track it.
			tot += (inst - da - db)
			cnt += 1.0
			
		# calculate average.
		return tot / cnt

			
	def _add_nodes(self):
		''' adds nodes from list.'''
		logging.debug("adding nodes.")
		
		# add nodes in one shot.
		self.add_nodes_from( list(self.node_list[:]['node_idx']) )


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
				subg.add_edge(edge[0],edge[1],weight=self[edge[0]][edge[1]]['gap_size'])
			
			# compute topological order.
			path = nx.topological_sort(subg)
			
			# validate path.
			pn = path[0]
			for i in range(1,len(path)):
				if subg.has_edge(pn, path[i]) == False:
					self.draw("failed2",subg)
					logging.error("Not a linear graph.")
					sys.exit()
				pn = path[i]
				
			
			# setup scaffold vars.
			scaf_idx = 1
			scaf_start = 1
			scaf_stop = self.node_list[path[0]]['ctg_width']
			
			# add first node.
			self.agp[agp_idx]['scaf_name'] = "scaf_%i" % idx
			self.agp[agp_idx]['scaf_start'] = scaf_start
			self.agp[agp_idx]['scaf_stop'] = scaf_stop
			self.agp[agp_idx]['scaf_idx'] = scaf_idx
			self.agp[agp_idx]['comp_type'] = "W"
			self.agp[agp_idx]['comp_name'] = self.node_list[path[0]]['ctg_name']
			self.agp[agp_idx]['comp_start'] = 1
			self.agp[agp_idx]['comp_stop'] = self.node_list[path[0]]['ctg_width']
			self.agp[agp_idx]['comp_orien'] = self.node_list[path[0]]['ctg_orien']
			self.agp[agp_idx]['comp_linkage'] = 0
			
			# increment pointers.
			agp_idx += 1
			pn = path[0]
			scaf_idx += 1
			scaf_start = scaf_stop + 1
			
			# fill AGP.
			for i in range(1,len(path)):

				# defualt gapsize to 7.
				gapsz = subg[pn][path[i]]['weight']
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
				scaf_stop = scaf_start + self.node_list[path[i]]['ctg_width']
				
				# add nodes.
				self.agp[agp_idx]['scaf_name'] = "scaf_%i" % idx
				self.agp[agp_idx]['scaf_start'] = scaf_start
				self.agp[agp_idx]['scaf_stop'] = scaf_stop
				self.agp[agp_idx]['scaf_idx'] = scaf_idx
				self.agp[agp_idx]['comp_type'] = "W"
				self.agp[agp_idx]['comp_name'] = self.node_list[path[i]]['ctg_name']
				self.agp[agp_idx]['comp_start'] = 1
				self.agp[agp_idx]['comp_stop'] = self.node_list[path[i]]['ctg_width']
				self.agp[agp_idx]['comp_orien'] = self.node_list[path[i]]['ctg_orien']
				
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
		
	def read_agp(self, fpath):
		
		# read in agp.
		fin = open(fpath, "rb")
		lines = fin.readlines()
		fin.close()

		# instantiate array.
		agp_edges = np.zeros(len(lines), dtype=agp_dt)
		
		# parse agp.
		idx = 0
		for line in lines:
			# tokenize.
			if line[0] == "#": continue
			tmp = line.strip().split()
			
			# get general tokenize.
			agp_edges[idx]['scaf_name'] = tmp[0]
			agp_edges[idx]['scaf_start'] = int(tmp[1])
			agp_edges[idx]['scaf_stop'] = int(tmp[2])
			agp_edges[idx]['scaf_idx'] = int(tmp[3])
			agp_edges[idx]['comp_type'] = tmp[4]
			
			# contig.
			if tmp[4] == "W":
				# get parts.
				agp_edges[idx]['comp_name'] = tmp[5]
				agp_edges[idx]['comp_start'] = int(tmp[6])
				agp_edges[idx]['comp_stop'] = int(tmp[7])
				if tmp[8] == "+":
					agp_edges[idx]['comp_orien'] = 0
				else:
					agp_edges[idx]['comp_orien'] = 1
				
			else:
				# save entry.
				agp_edges[idx]['comp_name'] = tmp[6]
				agp_edges[idx]['comp_start'] = 1
				agp_edges[idx]['comp_stop'] = int(tmp[5])
				if tmp[7] != "yes":
					agp_edges[idx]['comp_linkage'] = 0
				else:
					agp_edges[idx]['comp_linkage'] = 1
				

			# update index.
			idx += 1
		
		# shirnk array.
		agp_edges.resize(idx)
		
		# sort based on scaffold id, then part id.
		agp_edges.sort(order=['scaf_name','scaf_idx'])
		
		# create buffer for missing gaps.
		buffer = np.zeros(agp_edges.size / 2, dtype=agp_dt)
		buffer_idx = 0
		
		# validate and setup.
		scaf_name = ""
		scaf_idx = 0
		scaf_start = 0
		idx_bump = 0
		i = 0
		while i < agp_edges.size:
			
			# set variable.
			edge = agp_edges[i]
			
			# check if new entry.
			if edge['scaf_name'] != scaf_name:
				# validate new entry.
				if edge['scaf_idx'] != 1:
					logging.error("Scaffold not started at index 1.\n%s\n" % \
						(str(edge)))
					logging.error("bad1")
					sys.exit()

					
				if edge['comp_type'] != "W":
					logging.error("Scaffold not started with contig.\n%s\n" % \
						(str(edge)))
					logging.error("bad1")
					sys.exit()
				
				# setup new entry.
				scaf_name = edge['scaf_name']
				scaf_idx = edge['scaf_idx']
				scaf_start = edge['scaf_start']
				scaf_stop = edge['scaf_stop']
				gap_needed = True
				i += 1
				continue
			
			# check consistent idx.
			if edge['scaf_idx'] != scaf_idx + 1:
				print edge['scaf_idx'], scaf_idx + 1
				logging.error("Scaffold index out of order.\n%s\n" % \
					(str(edge)))
				logging.error("bad1")
				sys.exit()
				
			# check increasing start > stop.
			if edge['scaf_start'] >= edge['scaf_stop']:
				logging.error("Scaffold start stop not in order.\n%s\n" % \
					(str(edge)))
				logging.error("bad1")
				sys.exit()

			# check increasing start.
			if edge['scaf_start'] >= edge['scaf_stop']:
				logging.error("Component not in increasing order.\n%s\n" % \
					(str(edge)))
				logging.error("bad1")
				sys.exit()

			
			# check if we are looking for a gap.
			if gap_needed == True:
				
				# check if this is a gap.
				if edge['comp_type'] == "N":
									
					# load gap info.
					scaf_idx = edge['scaf_idx']
					scaf_start = edge['scaf_start']
					scaf_stop = edge['scaf_stop']
					gap_needed = False
					i += 1
					continue
					
				# file log.
				logging.warning("Missing gap.\n%s\n" % (str(edge)))
				#sys.exit()
				
			# handle continued component.
			scaf_idx = edge['scaf_idx']
			scaf_start = edge['scaf_start']
			scaf_stop = edge['scaf_stop']
			i += 1	
		
		
		# return array.
		return agp_edges		

	def draw(self, name, subg=False):
		return
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
