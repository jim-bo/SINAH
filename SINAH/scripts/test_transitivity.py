'''
looks for transitive reduction opportunities 
'''
# program imports.
from data_structs.nodes import load_nodes
from data_structs.nodes import create_lookup
from data_structs.edges import load_edges
from data_structs.bundles import load_bundles
from data_structs import types

# system imports.
import time
import random
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
input_nodes_file = os.path.abspath(sys.argv[1])
input_edges_file = os.path.abspath(sys.argv[2])

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

def get_triangles(bundles):
	''' gets triangles '''
	
	# build dict of neighbors.
	neibs = dict()
	bidxs = []
	for i in range(bundles.size):
		
		# simplify.
		idxa = bundles[i]['ctg_a_idx']
		idxb = bundles[i]['ctg_b_idx']
		
		
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
			
	return tris

########### script ################## 

# load hdf5 information.
logging.info("loading data arrays")
nodes = load_nodes(input_nodes_file)
edges = load_edges(input_edges_file)

# build edge lookup.
elookup = dict()
for i in range(edges.size):

	# get id
	idxa = edges[i]['ctg_a_idx']
	idxb = edges[i]['ctg_b_idx']

	elookup[(idxa, idxb)] = i
	elookup[(idxb, idxa)] = i
	
# get triangles.
tris = get_triangles(edges)
	
# search triangles for a small node.
for tri in tris:	
	for x in range(3):
		if nodes[tri[x]]['ctg_width'] < 1395 + (3 * 250):
				
			# look for triangle.
			i = tri[(x - 1) % 3]
			j = tri[x]
			k = tri[(x + 1) % 3]
			
			# check the edge from i to k can span 
			if (edges[elookup[(i,k)]]['implied_dist'] + (3 * 250)) > nodes[x]['ctg_width']:
				
				# simplify orientation.
				if edges[elookup[(i,j)]]['implied_state'] == 0 or edges[elookup[(i,j)]]['implied_state'] == 1:
					state_ij = 0
				else:
					state_ij = 1 
				if edges[elookup[(j,k)]]['implied_state'] == 0 or edges[elookup[(j,k)]]['implied_state'] == 1:
					state_jk = 0
				else:
					state_jk = 1
				if edges[elookup[(i,k)]]['implied_state'] == 0 or edges[elookup[(i,k)]]['implied_state'] == 1:
					state_ik = 0
				else:
					state_ik = 1
				
				# check the orientation.
				if state_ij == state_jk and state_jk == state_ik:
					print i,j,k
					break
	
	
	
