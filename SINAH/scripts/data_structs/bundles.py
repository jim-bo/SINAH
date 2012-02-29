'''
functions to create edge and edge weights in python
and saving the files to HDF5 tile.
'''
# program imports.
from data_structs.types import bundle_dt

# system imports.
import numpy as np
import logging
import h5py
import sys

########### functions ##################

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

def load_bundles(file_path):
	''' loads edges from h5py file.'''
	'''
	# open the file read only.
	f = h5py.File(file_path, 'r')
	
	# load data into memory.
	data = f['bundles'][:]
		
	# close file.
	f.close()
	'''
	
	# load files into memory.
	fin = open(file_path, "rb")
	lines = fin.readlines()
	fin.close()
	
	# create numpy array.
	data = np.zeros(len(lines), dtype=bundle_dt)
	
	# populate it.
	i = 0
	for line in lines:
		line = line.split("\t")
		data[i]['ctg_a_idx'] = int(line[0])
		data[i]['ctg_b_idx'] = int(line[1])
		data[i]['WT_A'] = float(line[2])
		data[i]['WT_B'] = float(line[3])
		data[i]['WT_C'] = float(line[4])
		data[i]['WT_D'] = float(line[5])
		i += 1
	
	
	# return data.
	return data

def save_bundles(bundles, file_path):
	''' save edge array ot hd5f file.'''
	
	# open the file.
	f = h5py.File(file_path, 'w')
	
	# save the dataset.
	f.create_dataset('bundles', data=bundles)
	
	# close the file.
	f.close()
	

