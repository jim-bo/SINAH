'''
functions to create edge and edge weights in python
and saving the files to HDF5 tile.
'''
# program imports.
from data_structs.types import node_dt
from data_structs.types import edge_dt

# system imports.
import subprocess
import numpy as np
import logging
import h5py
import sys
import mmap
import gzip

########### functions ##################

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
		sys.exit(1)
	if d2 < 0:
		print "bad dist: d2"
		print d2
		print state
		print aid, bid
		sys.exit(1)
		
	# finalize dist.
	dist = edge['insert_size'] - d1 - d2
		
	# return dist.
	return dist
		
		
	

def determine_class(aid, bid, ao, bo):
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

def read_sam_skim(sam_file, id_key, cutoff=-1):
	''' preprocess sam file to find pairs.'''
	
	# read same file.
	if sam_file.count(".gz") > 0:
		fin = gzip.open(sam_file)
	else:
		fin = open(sam_file, "rb")
	
	# parse.
	debug = 0
	offset = 0
	res = {}
	for line in fin:
		# skips.
		if line[0] == "@": 
			offset += len(line)
			continue
			
		tmp = line.split("\t")
		
		if len(tmp) <= 2:
			offset += len(line)
			continue
		
		if tmp[2] == "*": 
			offset += len(line)
			continue
			
		# tokenize.
		id = tmp[0].replace(id_key,"")
		
		# check if not fresh.
		if id in res:
			res[id] = -1
			offset += len(line)
			continue
			
		# insert.
		res[id] = offset
		
		# increment offset.
		offset += len(line)
		
		if cutoff != -1:
			debug += 1
			if debug > cutoff:
				break
		
	# return result.
	fin.close()
	return res

def load_edges(file_path):
	''' loads edges from h5py file.'''
	
	'''
	# open the file read only.
	f = h5py.File(file_path, 'r')
	
	# load data into memory.
	data = f['edges'][:]
		
	# close file.
	f.close()
	'''
	
	# load files into memory.
	fin = open(file_path, "rb")
	lines = fin.readlines()
	fin.close()
	
	# create numpy array.
	data = np.zeros(len(lines), dtype=edge_dt)
	
	# populate it.
	i = 0
	for line in lines:
		line = line.split("\t")
		data[i]['ctg_a_idx'] = int(line[0])
		data[i]['ctg_b_idx'] = int(line[1])
		data[i]['read_a_left_pos'] = int(line[2])
		data[i]['read_a_right_pos'] = int(line[3])
		data[i]['read_b_left_pos'] = int(line[4])
		data[i]['read_b_right_pos'] = int(line[5])
		data[i]['read_a_orien'] = int(line[6])
		data[i]['read_b_orien'] = int(line[7])
		data[i]['insert_size'] = int(line[8])
		data[i]['implied_state'] = int(line[9])
		data[i]['implied_dist'] = int(line[10])
		data[i]['std_dev'] = int(line[11])
		data[i]['invalid'] = int(line[12])
		data[i]['used'] = int(line[13])
		i += 1
	
	# return data.
	return data

def save_edges(edges, file_path):
	''' save edge array ot hd5f file.'''
	
	logging.error("not supported")
	sys.exit(1)	
	
	subprocess.call(["rm","-rf",file_path])
	
	# open the file.
	f = h5py.File(file_path, 'w')
	
	# save the dataset.
	f.create_dataset('edges', data=edges)
	
	# close the file.
	f.close()
	
