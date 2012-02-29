'''
functions to create edge and edge weights in python
and saving the files to HDF5 tile.
'''
# program imports.
from data_structs.types import bundle_dt

# system imports.
import biobuffers
import numpy as np
import logging
import h5py
import sys

########### functions ##################

def create_bundles(nodes, edges, weights):
	''' creates bundle given marked edges.'''
	
	# determine number of bundles.
	keys = {}
	
	# loop over edges.
	for i in range(edges.size):
		
		# skip invalid.
		if edges[i]['invalid'] == True:
			continue
		
		# make key.
		key = tuple(sorted([edges[i]['ctg_a_idx'], edges[i]['ctg_b_idx']]))
		
		# add index to set.
		if key not in keys:
			keys[key] = []
		keys[key].append(i)	
		
	# allocate array for bundles.
	sz = len(keys)
	bundles = np.zeros(sz, dtype=bundle_dt)
	
	# copy data into array.
	i = 0
	for key in keys:
		
		# build probabilities.
		wt_A = 0.0
		wt_B = 0.0
		wt_C = 0.0
		wt_D = 0.0
		
		# loop over indices.
		state_lists = ([], [], [], [])
		for j in keys[key]:
			
			# calculate pk.
			pk = weights[j]['a_overlap'] * weights[j]['b_overlap']
			
			# get vars.
			aid = edges[j]['ctg_a_idx']
			bid = edges[j]['ctg_b_idx']
			ao = edges[j]['read_a_orien']
			bo = edges[j]['read_b_orien']
			
			# determine class. 
			type = determine_class(aid, bid, ao, bo)
			
			# add weight based on class.
			if type == 0:
				wt_A += pk
				state_lists[0].append(j)
			elif type == 1:
				wt_B += pk
				state_lists[1].append(j)
			elif type == 2:
				wt_C += pk
				state_lists[2].append(j)
			elif type == 3:
				wt_D += pk
				state_lists[3].append(j)
		
		# identify largest window.
		'''
		state_wts = [[], [], [], []]
		for q in range(4):
					
			# build array of start positions.
			#tmp1 = edges[state_lists[i]]['read_a_left_pos']
			tmp = edges[state_lists[q]]['implied_dist']
			
			# sanity.
			if len(tmp) < 2:
				continue
				
			# sort the array indirectly.
			tmp_keys = np.argsort(tmp)
			
			# find max window.
			best_win = []
			best_score = -1
			ltmp = len(tmp)
			for j in range(ltmp):
				
				# loop until 2 * 3 * std_dev is scanned.
				cur_win = [tmp[tmp_keys[j]]]
				for k in range(j+1,ltmp):
					
					print tmp[tmp_keys[j]]
					
					# add score to window.
					if (tmp[tmp_keys[k]] - tmp[tmp_keys[j]]) <= (2 * 3 * edges[tmp[tmp_keys[j]]]['std_dev']):
						cur_win.append(tmp[tmp_keys[k]])
					else:
						break
					
				# score it.
				score = 1.0
				for k in range(len(cur_win)):
					score += cur_win[k]
				
				# check if replace.
				if score > best_score:
					best_score = score
					best_win = cur_win
					
			# update state weight.
			state_wts[q] = best_win

		'''		
		# save ids to bundles.
		bundles[i]['ctg_a_idx'] = key[0]
		bundles[i]['ctg_b_idx'] = key[1]
				
		# save counts to bundles.
		#'''
		bundles[i]['WT_A'] = wt_A
		bundles[i]['WT_B'] = wt_B
		bundles[i]['WT_C'] = wt_C
		bundles[i]['WT_D'] = wt_D
		#'''
		'''
		bundles[i]['WT_A'] = len(state_wts[0])
		bundles[i]['WT_B'] = len(state_wts[1])
		bundles[i]['WT_C'] = len(state_wts[2])
		bundles[i]['WT_D'] = len(state_wts[3])
		'''	
		'''
		prob_wt = [1, 1, 1, 1]
		for j in range(4):
			for k in range(len(state_wts[j])):
				prob_wt[j] += state_wts[j][k]
		bundles[i]['WT_A'] = prob_wt[0]
		bundles[i]['WT_B'] = prob_wt[1]
		bundles[i]['WT_C'] = prob_wt[2]
		bundles[i]['WT_D'] = prob_wt[3]
		'''
		#print wt_A, wt_B, wt_C, wt_D, state_wts
		
		# increment bundle idx.
		i += 1	
			
		
	# return bundles.
	return bundles

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
	

