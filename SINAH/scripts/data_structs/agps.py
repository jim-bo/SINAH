'''
functions to create edge and edge weights in python
and saving the files to HDF5 tile.
'''
# program imports.
from data_structs.types import agp_dt

# system imports.
import subprocess
import numpy as np
import logging
import sys

########### functions ##################

def load_agps(fpath):
	''' read agp file into array.'''
	
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
		if tmp[0] == "RUNTIME:": continue
		
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
	
	return agp_edges


def save_agps(agp_file, agp):
	''' saves agp to disk.'''
	
	# write to file.
	fout = open(agp_file, "w")
	
	# write each entry.
	z = len(agp_dt.names)
	for i in range(agp.size):
		
		# format result.
		tmp = agp[i]
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
