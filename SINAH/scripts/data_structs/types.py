'''
Created on Apr 8, 2011

data types used in the graph.

@author: eljimbo
'''

import numpy as np

node_dt = np.dtype([\
	('node_idx', np.long),\
	('ctg_name', 'S255'),\
	('ctg_width', np.long),\
	('ctg_orien', np.long),\
	('ctg_order', np.long),\
	('invalid', np.long)\
])

edge_dt = np.dtype([\
	('ctg_a_idx', np.long),\
	('ctg_b_idx', np.long),\
	('read_a_left_pos', np.long),\
	('read_a_right_pos', np.long),\
	('read_b_left_pos', np.long),\
	('read_b_right_pos', np.long),\
	('read_a_orien', np.long),\
	('read_b_orien', np.long),\
	('insert_size', np.long),\
	('implied_state', np.long),\
	('implied_dist', np.long),\
	('std_dev', np.long),\
	('invalid', np.long),\
	('used', np.long)\
])

bundle_dt = np.dtype([\
	('ctg_a_idx', np.long),\
	('ctg_b_idx', np.long),\
	('WT_A', np.double),\
	('WT_B', np.double),\
	('WT_C', np.double),\
	('WT_D', np.double)
])

directed_dt = np.dtype([\
	('a_idx', np.long),\
	('b_idx', np.long),\
	('a_d', np.long),\
	('b_d', np.long),\
	('gap_size', np.long),\
	('insert_size', np.long)
])

agp_dt = np.dtype([\
	('scaf_name', 'S255'),\
	('scaf_start', np.long),\
	('scaf_stop', np.long),\
	('scaf_idx', np.long),\
	('comp_type', 'S50'),\
	('comp_name', 'S255'),\
	('comp_start', np.long),\
	('comp_stop', np.long),\
	('comp_orien', np.long),\
	('comp_linkage', np.long),\
])

weight_dt = np.dtype([\
	('a_idx', np.long),\
	('b_idx', np.long),\
	('a_coverage', np.double),\
	('b_coverage', np.double),\
	('a_overlap', np.double),\
	('b_overlap', np.double),\
])
