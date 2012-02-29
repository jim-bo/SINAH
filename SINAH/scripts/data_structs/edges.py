'''
functions to create edge and edge weights in python
and saving the files to HDF5 tile.
'''
# program imports.
from data_structs.types import node_dt
from data_structs.types import edge_dt

# system imports.
import subprocess
import biobuffers
import numpy as np
import logging
import h5py
import sys
import mmap
import gzip

########### functions ##################

def create_edges_new(nodes, node_lookup, sam_file1, sam_file2, insert_size, std_dev, cutoff=-1):
	''' creates edges given two SAM files.'''
	
	# find number of edges.
	fin1 = open(sam_file1, "rb")
	num_edges = 0
	for line in fin1.readlines():
		num_edges += 1
	fin1.close()
	
	# create array for edges.
	logging.debug("allocating edge array for %i putative edges." % num_edges)
	edges = np.zeros(num_edges, dtype=edge_dt)
	
	# read in valid SAM entries.
	fin1 = open(sam_file1, "rb")
	fin2 = open(sam_file2, "rb")
	j = 0
	same_ctg = 0
	logging.debug("traversing array")
	for i in range(edges.size):
		
		if i % 1000000 == 0:
			logging.debug("... %i of %i" % (i, num_edges))
		
		# read lines.
		line1 = fin1.readline()
		line2 = fin2.readline()
		
		# tokenize lines.
		tmp1 = line1.split("\t")
		tmp2 = line2.split("\t")

		# check dictionary for contig id.
		if tmp1[2] not in node_lookup:
			continue
			logging.error("contig missing in node file: %s" % tmp1[2])
			sys.exit()

		if tmp2[2] not in node_lookup:
			continue
			logging.error("contig missing in node file: %s" % tmp2[2])
			sys.exit()
						
		id1 = node_lookup[tmp1[2]]
		id2 = node_lookup[tmp2[2]]

		# note if they map to same contig.
		if id1 == id2:
			same_ctg += 1
			continue

		# parse orientation.
		if tmp1[1] == "0":  o1 = 0
		else:			   o1 = 1
		
		if tmp2[1] == "0":  o2 = 0
		else:			   o2 = 1	

		# add info to entry.
		edges[j]['ctg_a_idx'] = id1
		edges[j]['ctg_b_idx'] = id2
		edges[j]['read_a_left_pos'] = int(tmp1[3])
		edges[j]['read_a_right_pos'] = int(tmp1[3]) + len(tmp1[9])
		edges[j]['read_b_left_pos'] = int(tmp2[3]) 
		edges[j]['read_b_right_pos'] = int(tmp2[3]) + len(tmp2[9])
		edges[j]['read_a_orien'] = o1
		edges[j]['read_b_orien'] = o2
		edges[j]['insert_size'] = insert_size
		edges[j]['std_dev'] = std_dev
		edges[j]['invalid'] = False

		# determine state.
		state = determine_class(id1, id2, o1, o2)
		
		# determine distance.
		dist = determine_dist(id1, id2, state, edges[j], nodes[id1], nodes[id2])

		# save state and dist to edge.
		edges[j]['implied_state'] = state
		edges[j]['implied_dist'] = dist
		
		# increment j.
		j+= 1
	 
	# close read files.
	fin1.close()
	fin2.close()
	
	# resize array to appropriate size.
	logging.debug("found %i same ctg edges" % same_ctg)
	logging.debug("found %i valid edges" % j)
	edges.resize(j)
	
	# return edges
	return edges

def create_edges(nodes, node_lookup, sam_file1, sam_file2, id1, id2, insert_size, std_dev, cutoff=-1):
	''' creates edges given two SAM files.'''
	
	# read first file.
	logging.debug("reading first")
	mapped_1 = read_sam_skim(sam_file1, id1, cutoff=cutoff)
	
	# read second file.
	logging.debug("reading second")
	mapped_2 = read_sam_skim(sam_file2, id2, cutoff=cutoff)
	
	# identify list of usable pairs.
	usable = []
	for id in mapped_1:
		
		# check if unique.
		if mapped_1[id] == -1: continue
		
		# check if in other.
		if id not in mapped_2: continue
		
		# check if unique in other.
		if mapped_2[id] == -1: continue
		
		# add pair info.
		usable.append(id)
		
	# save size.
	sz = len(usable)
		
		
	# arrays to track order.
	logging.debug("populating order array")
	tmp_dt = np.dtype( [ ('id', 'S100'), ('idx1', np.int), ('idx2', np.int) ] )
	orders = np.zeros(sz, dtype=tmp_dt)
	
	# populate array with entries.
	i = 0
	for id in usable:
		
		if i % 1000000 == 0:
			logging.debug("... %i of %i" % (i, sz))
		
		# get idx.
		idx1 = mapped_1[id]
		idx2 = mapped_2[id]
		
		# add to line.
		orders[i]['id'] = id
		orders[i]['idx1'] = idx1
		orders[i]['idx2'] = idx2
		i += 1
		
	# sort by idx1 then idx2.
	logging.debug("sorting order array")
	orders.sort(order=['idx1','idx2'])
		
	# create array for edges.
	logging.debug("allocating edge array for %i putative edges." % sz)
	edges = np.zeros(sz, dtype=edge_dt)
	
	# count faults
	same_ctg = 0
	
	# read in valid SAM entries.
	fin1 = open(sam_file1, "rb")
	fin2 = open(sam_file2, "rb")
	j = 0
	logging.debug("traversing array")
	for i in range(orders.size):
		
		if i % 1000000 == 0:
			logging.debug("... %i of %i" % (i, sz))
		
		# seek to lines.
		fin1.seek(orders[i]['idx1'])
		fin2.seek(orders[i]['idx2'])
		line1 = fin1.readline()
		line2 = fin2.readline()
		
		
		# tokenize lines.
		tmp1 = line1.split("\t")
		tmp2 = line2.split("\t")


		# check dictionary for contig id.
		if tmp1[2] not in node_lookup:
			continue
			logging.error("contig missing in node file: %s" % tmp1[2])
			sys.exit()

		if tmp2[2] not in node_lookup:
			continue
			logging.error("contig missing in node file: %s" % tmp2[2])
			sys.exit()
						
		id1 = node_lookup[tmp1[2]]
		id2 = node_lookup[tmp2[2]]

		# note if they map to same contig.
		if id1 == id2:
			same_ctg += 1
			continue

		# parse orientation.
		if tmp1[1] == "0":  o1 = 0
		else:			   o1 = 1
		
		if tmp2[1] == "0":  o2 = 0
		else:			   o2 = 1	

		# add info to entry.
		edges[j]['ctg_a_idx'] = id1
		edges[j]['ctg_b_idx'] = id2
		edges[j]['read_a_left_pos'] = int(tmp1[3])
		edges[j]['read_a_right_pos'] = int(tmp1[3]) + len(tmp1[9])
		edges[j]['read_b_left_pos'] = int(tmp2[3]) 
		edges[j]['read_b_right_pos'] = int(tmp2[3]) + len(tmp2[9])
		edges[j]['read_a_orien'] = o1
		edges[j]['read_b_orien'] = o2
		edges[j]['insert_size'] = insert_size
		edges[j]['std_dev'] = std_dev
		edges[j]['invalid'] = False

		# determine state.
		state = determine_class(id1, id2, o1, o2)
		
		# determine distance.
		dist = determine_dist(id1, id2, state, edges[j], nodes[id1], nodes[id2])

		# save state and dist to edge.
		edges[j]['implied_state'] = state
		edges[j]['implied_dist'] = dist
		
		# increment j.
		j+= 1
	 
	# close read files.
	fin1.close()
	fin2.close()
	
	# resize array to appropriate size.
	logging.debug("found %i same ctg edges" % same_ctg)
	logging.debug("found %i valid edges" % j)
	edges.resize(j)
	
	# return edges
	return edges

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
	
