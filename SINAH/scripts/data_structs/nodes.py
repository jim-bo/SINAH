'''
functions to create nodes in python
and saving the files to HDF5 tile.
'''
# program imports.
from data_structs.types import node_dt

# system imports.
import biobuffers
import numpy as np
import logging
import h5py
import sys

############## functions.##################

def create_lookup(nodes):
	''' creates lookup dictionary '''
	# loop over names.
	lookup = {}
	
	for node in nodes:
		lookup[node['ctg_name']] = node['node_idx']
		
	return lookup
	

def create_nodes(fasta_file):
	''' creates node array given fasta file.'''
	
	# loop over each seq.
	logging.debug("retrieving nodes information.")
	info = []
	idx = 0
	for entry in biobuffers.fasta_buffer(fasta_file):
		
		# save info to list.
		info.append((idx, entry['name'], len(entry['seq'])))
		idx += 1
	
	
	# create array for nodes.
	sz = len(info)
	logging.debug("creating array for %i nodes." % sz)
	nodes = np.zeros(sz, dtype=node_dt)
	
	# save list to structure.
	logging.debug("copying %i nodes." % sz)
	for i in range(nodes.size):
		# debug.
		if i % 100000 == 0:
			logging.debug("... %i" % i)
		
		# set info.
		nodes[i]['node_idx'] = info[i][0]
		nodes[i]['ctg_name'] = info[i][1]
		nodes[i]['ctg_width'] = info[i][2]
		nodes[i]['ctg_orien'] = -1
		nodes[i]['ctg_order'] = -1
		
				
	# return list.
	return nodes

def load_nodes(file_path):
	''' loads nodes from h5py file.'''
	
	'''
	# open the file read only.
	f = h5py.File(file_path, 'r')
	
	# load data into memory.
	data = f['nodes'][:]
	'''
	
	# load files into memory.
	fin = open(file_path, "rb")
	lines = fin.readlines()
	fin.close()
	
	# create numpy array.
	data = np.zeros(len(lines), dtype=node_dt)
	
	# populate it.
	i = 0
	for line in lines:
		line = line.split("\t")
		data[i]['node_idx'] = int(line[0])
		data[i]['ctg_name'] = line[1]
		data[i]['ctg_width'] = int(line[2])
		data[i]['ctg_orien'] = int(line[3])
		data[i]['ctg_order'] = int(line[4])
		data[i]['invalid'] = int(line[5])
		i += 1
		
	# return data.
	return data
		

def save_nodes(nodes, file_path):
	''' save node array ot hd5f file.'''
	logging.error("not supported")
	sys.exit(1)
	
	'''
	# open the file.
	f = h5py.File(file_path, 'w')
	
	# save the dataset.
	f.create_dataset('nodes', data=nodes)
	
	# close the file.
	f.close()
	'''
