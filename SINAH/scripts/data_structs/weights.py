'''
functions to create edge and edge weights in python
and saving the files to HDF5 tile.
'''
# program imports.
from data_structs.types import node_dt
from data_structs.types import edge_dt
from data_structs.types import weight_dt

# system imports.
from biobuffers import gff_buffer
from bitarray import bitarray

import numpy as np
import logging
import h5py
import sys

########### functions ##################

def create_weights(nodes, node_lookup, edges, gff_file):
	''' creates weights from gff file.'''

	# allocate weight structure.
	logging.debug("allocating weight array")
	weights = np.zeros(edges.size, dtype=weight_dt)
	
	# create bit dictionary.
	logging.debug("creating bit arrays.")
	bit_dict = create_dict(nodes)

	# load annotation from gff.
	logging.debug("annotating bit arrays.")
	annotate_genome(gff_file, bit_dict)

	# loop over edges.
	logging.debug("processing edges")
	for i in range(edges.size):

		if i % 100000 == 0:
			logging.debug("... %i " % i)

		# get contig names.
		edge = edges[i]
		ctg1_name = nodes[edge['ctg_a_idx']]['ctg_name']
		ctg2_name = nodes[edge['ctg_b_idx']]['ctg_name']

		# set idx vars.
		weights[i]['a_idx'] = edge['ctg_a_idx']
		weights[i]['b_idx'] = edge['ctg_b_idx']
			
		# add coverage to array, will be turned to prob later.
		weights[i]['a_coverage'] = 1.0
		weights[i]['b_coverage'] = 1.0
		
		# get overlap.
		r1_olap = bit_dict[ctg1_name][edge['read_a_left_pos']:edge['read_a_right_pos']].count(1)
		r2_olap = bit_dict[ctg2_name][edge['read_b_left_pos']:edge['read_b_right_pos']].count(1)
		
		# add overlap probability.
		r1_sz = edge['read_a_right_pos'] - edge['read_a_left_pos']
		r2_sz = edge['read_b_right_pos'] - edge['read_b_left_pos']
		
		weights[i]['a_overlap'] = 1.0 - ( float(r1_olap) / float(r1_sz) )
		weights[i]['b_overlap'] = 1.0 - ( float(r2_olap) / float(r2_sz) )
	
	# return weights.
	return weights
	
def create_dict(nodes):
	''' creates bit array dictionary '''
	
	# create dict.
	bit_dict = {}
	
	# loop over each contig.
	for node in nodes:
		
		# create array.
		tmp = bitarray(node['ctg_width'])
		
		# initialize to zero.
		tmp[:] = False
		
		# save to dictionary.
		bit_dict[node['ctg_name']] = tmp
		
	# return result.
	return bit_dict

def annotate_genome(gff_file, bit_dict):
	''' marks filtered regions with 1'''
		
	# loop over entries.
	cnt = 0
	for entry in gff_buffer(gff_file):
			
		# skip bad entries.
		if entry['seqname'] not in bit_dict:
			continue
			
		# mark as repetative.
		bit_dict[entry['seqname']][entry['start']:entry['end']] = True


def load_weights(file_path):
	''' loads weights from h5py file.'''
	
	# open the file read only.
	f = h5py.File(file_path, 'r')
	
	# load data into memory.
	data = f['weights'][:]
		
	# return data.
	return data

def save_weights(weights, file_path):
	''' save edge array ot hd5f file.'''
	
	# open the file.
	f = h5py.File(file_path, 'w')
	
	# save the dataset.
	f.create_dataset('weights', data=weights)
	
	# close the file.
	f.close()
	
