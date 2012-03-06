'''
breaks large 2-components further.
'''
# program imports.
from data_structs.nodes import load_nodes
from data_structs.nodes import create_lookup
from data_structs.bundles import load_bundles
from data_structs import types

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

print "AAAAAAAAAAAA"
sys.exit()

# logging.
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )
#logging.basicConfig(level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )

# parameters.
input_nodes_file = os.path.abspath(sys.argv[1])
input_bundles_file = os.path.abspath(sys.argv[2])
graph_dir = os.path.abspath(sys.argv[3])
cutoff = int(sys.argv[4])

decomp_0_file = "%s/decomp_zero.txt" % graph_dir
decomp_1_file = "%s/decomp_one.txt" % graph_dir
decomp_2_file = "%s/decomp_two.txt" % graph_dir


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

# load hdf5 information.
logging.info("loading data arrays")
nodes = load_nodes(input_nodes_file)
bundles = load_bundles(input_bundles_file)

# load decomposition.
DG = load_decomposition(decomp_0_file, decomp_1_file, decomp_2_file)

# build the neighbor list.
neighbors = dict()
for i in range(bundles.size):
	
	ida = bundles[i]['ctg_a_idx']
	idb = bundles[i]['ctg_b_idx']
	
	if ida not in neighbors:
		neighbors[ida] = list()
	neighbors[ida].append(idb)

	if idb not in neighbors:
		neighbors[idb] = list()
	neighbors[idb].append(ida)
	

# find too large components.
for n0 in DG.nodes():
	
	# skip.
	if len(DG.node[n0]['set']) < cutoff: continue
	
	# next level.
	for n1 in DG.node[n0]['graph'].nodes():

		# skip.
		if len(DG.node[n0]['graph'].node[n1]['set']) < cutoff: continue
		
		# next level
		for n2 in DG.node[n0]['graph'].node[n1]['graph'].nodes():
			
			# skip.
			if len(DG.node[n0]['graph'].node[n1]['graph'].node[n2]['set']) < cutoff: continue
			
			# simplify.
			s = DG.node[n0]['graph'].node[n1]['graph'].node[n2]['set']
			
			# map all nodes.
			amap = dict()
			bmap = dict()
			idx = 1
			for n in s:
				amap[n] = idx
				bmap[idx] = n
				idx += 1
			
			# build induced graph.
			g = nx.Graph()
			for i in range(bundles.size):
				idx1 = bundles[i]['ctg_a_idx']
				idx2 = bundles[i]['ctg_b_idx']
				if idx1 in s and idx2 in s:
					g.add_edge(amap[idx1], amap[idx2])
			
			# open file.
			fout = open("metis.in", "wb")
			
			# write header.
			fout.write("%i %i\n" % (g.number_of_nodes(), g.number_of_edges()))
				
			# write out info.
			node_cnt = len(g.nodes())
			for p in sorted(g.nodes()):
				for q in g.neighbors(p):
					fout.write("%i " % q)
				fout.write("\n")

			# close file.
			fout.close()
				
			# sanity check.		
			if nx.is_connected(g) == False:
				print "TROLOLOLOLO"
				sys.exit()
				
						
			# call metis.
			way_cut = node_cnt / cutoff
			
			print node_cnt, cutoff, way_cut
			
			sys.exit()
			subprocess.call(["gpmetis", "-objtype=vol", "metis.in", str(way_cut)])
			
			# open result.
			fin = open("metis.in.part.%i" % way_cut, "rb")
			
			# load partition.
			psets = dict()
			npart = dict()
			idx = 1
			for line in fin:
				# tokenize.
				n = bmap[idx]
				p = int(line.strip())
				
				# build the set.
				if p not in psets:
					psets[p] = set()
				psets[p].add(n)
				
				# build node lookup.
				npart[n] = p 
				
				# get next node.
				idx += 1

				
			# close file.
			fin.close()
				
			# fill in neighbors.
			for n in npart:
				for m in neighbors[n]:
					psets[npart[n]].add(m)
				
				
			# determine cut.
			for p1 in psets:
				for p2 in psets:
					if p1 == p2: continue
					
					# compare sets.
					intr = psets[p1].intersection(psets[p2])
					z = len(intr)
					
					# print cut.
					if z > 0:
						print p1, p2, len(intr)
				
			

				
			

			sys.exit()
				
