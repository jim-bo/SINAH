'''
tracks SPQR solution.

don't call this directly.
'''
# program imports.
from data_structs.agps import load_agps

# system imports.
import numpy as np
import time
import logging
import sys
import os

# logging.
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )
#logging.basicConfig(level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )

sol_node_dt = np.dtype([\
	('idx', np.int),\
	('orien', np.int)\
])

sol_bundle_dt = np.dtype([\
	('idxa', np.int),\
	('idxb', np.int),\
	('S', np.int),\
	('A', np.int),\
	('B', np.int),\
	('C', np.int),\
	('D', np.int),\
	('X', np.int),\
	('Y', np.int)\
])

class PartialSolution(object):
	''' paritial solution to ILP dictionary '''
	
	def __init__(self, szn, szb):
		''' create array with dictionary access '''
		
		# save size.
		self._nsz = szn
		self._bsz = szb
		
		# allocate arrays.
		self._node_tmp = np.zeros(self._nsz, dtype=sol_node_dt)
		self._bundle_tmp = np.zeros(self._bsz, dtype=sol_bundle_dt)
		
		# bootstrap index vars.
		self._node_idx = 0
		self._bundle_idx = 0

	def get_nsize(self):
		return self._node_tmp.size

	def get_bsize(self):
		return self._bundle_tmp.size
		
	def get_nodes(self):
		return self._node_tmp
		
	def get_bundles(self):
		return self._bundle_tmp


	def set_node(self, idx, val):
		''' set orientation '''
		
		# sanity check.
		if self._node_idx >= self._node_tmp.size:
			logging.error("bad sync between solution object and ILP: nodes: %i %i" % (self._node_tmp.size, self._node_idx))
			sys.exit(1)
		
		# set the value.
		self._node_tmp[self._node_idx]['idx'] = idx
		self._node_tmp[self._node_idx]['orien'] = val
		
		# update tracker.
		self._node_idx += 1
		
	def set_bundle(self, idxa, idxb, S, A, B, C, D, X, Y):
		''' set other vars '''
		
		# sanity check.
		if self._bundle_idx >= self._bundle_tmp.size:
			logging.error("bad sync between solution object and ILP: bundle")
			sys.exit(1)		
		
		# set the values.
		self._bundle_tmp[self._bundle_idx]['idxa'] = idxa
		self._bundle_tmp[self._bundle_idx]['idxb'] = idxb
		self._bundle_tmp[self._bundle_idx]['S'] = S
		self._bundle_tmp[self._bundle_idx]['A'] = A
		self._bundle_tmp[self._bundle_idx]['B'] = B
		self._bundle_tmp[self._bundle_idx]['C'] = C
		self._bundle_tmp[self._bundle_idx]['D'] = D
		self._bundle_tmp[self._bundle_idx]['X'] = X
		self._bundle_tmp[self._bundle_idx]['Y'] = Y
		
		
		# update the tracker.
		self._bundle_idx += 1

class SpqrSolution(object):
	''' tracks solution to scaffolding ILP '''
	
	def __init__(self, nodes, bundles):
		''' sets up solution arrays '''
		
		# grab size information.
		self._node_size = nodes.size
		self._bundle_size = bundles.size
		
		# allocate solution arrays.
		self._sol_nodes = np.zeros(self._node_size, dtype=sol_node_dt)
		self._sol_bundles = np.zeros(self._bundle_size, dtype=sol_bundle_dt)
		
		# default solutions to -1.
		self._sol_nodes[:]['idx'] = -1
		self._sol_nodes[:]['orien'] = -1
		self._sol_bundles[:]['idxa'] = -1
		self._sol_bundles[:]['idxb'] = -1
		self._sol_bundles[:]['S'] = -1
		self._sol_bundles[:]['A'] = -1
		self._sol_bundles[:]['B'] = -1
		self._sol_bundles[:]['C'] = -1
		self._sol_bundles[:]['D'] = -1
		
		# setup tracking vars.
		self._sol_added = 0
		self._nodes_added = set()
		
		# build node lookup dict.
		self._node_lookup = dict()
		for i in range(self._node_size):
			self._node_lookup[nodes[i]['ctg_name']] = nodes[i]['node_idx']
			
		# build bundle lookup dict.
		self._bundle_lookup = dict()
		for i in range(self._bundle_size):
			self._bundle_lookup[self._bundle_key(bundles[i]['ctg_a_idx'], bundles[i]['ctg_b_idx'])] = i

	def get_bslice(self, isect):
		''' return slice of nodes given a set '''
		
		# find edges with both.
		alist = []
		blist = []
		for i in range(self._sol_bundles.size):
			
			# check for both index.
			if self._sol_bundles[i]['idxa'] in isect or self._sol_bundles[i]['idxb'] in isect:
				alist.append(i)
				
		# return slice.
		return self._sol_bundles[alist]
		
	def get_orien(self, i):
		''' return orientation given an index '''
		
		# sanity check.
		if self._sol_nodes[i]['orien'] == -1:
			logging.error("getting orientation of unsolved node")
			
		return self._sol_nodes[i]['orien']

	def get_nslice(self, isect):
		''' return slice of nodes given a set '''
		return self._sol_nodes[list(isect)]

	def get_intersection(self, active):
		''' determines intersection '''
		return self._nodes_added.intersection(active)

	def get_bundles(self):
		''' return solution '''
		return self._sol_bundles

	def get_nodes(self):
		''' return solution '''
		return self._sol_nodes
		

	def _bundle_key(self, a, b):
		''' makes a bundle key '''
		
		# sort by largets first.
		if a > b:
			return (a, b)
		else:
			return (b, a)
		
	def _just_nodes(self, pnodes):
		''' adds nodes with no care in the world '''
	
		# add the node info.
		for i in range(pnodes.size):
			
			# get index.
			idx = pnodes[i]['idx']
			
			# save info.
			self._sol_nodes[idx]['idx'] = idx
			self._sol_nodes[idx]['orien'] = pnodes[i]['orien']
	
	def _just_bundles(self, pbundles):
		''' adds bundles with no care in the world '''

		# add the bundle info.
		for i in range(pbundles.size):
			
			# get index.
			idxa = pbundles[i]['idxa']
			idxb = pbundles[i]['idxb']
			idxbun = self._bundle_lookup[self._bundle_key(idxa, idxb)]
			
			# add bundles.
			self._sol_bundles[idxbun]['idxa'] = pbundles[i]['idxa']
			self._sol_bundles[idxbun]['idxb'] = pbundles[i]['idxb']
			
			self._sol_bundles[idxbun]['A'] = pbundles[i]['A']
			self._sol_bundles[idxbun]['B'] = pbundles[i]['B']
			self._sol_bundles[idxbun]['C'] = pbundles[i]['C']
			self._sol_bundles[idxbun]['D'] = pbundles[i]['D']
			self._sol_bundles[idxbun]['S'] = pbundles[i]['S']
			self._sol_bundles[idxbun]['X'] = pbundles[i]['X']
			self._sol_bundles[idxbun]['Y'] = pbundles[i]['Y']
	
	def apply_tri_partial(self, partial):
		''' may require flipping '''
		
		# try normal.
		status = self.apply_partial(partial, tri=True)
		if status == True:
			return
			
		# flip everything.
		pnodes = partial.get_nodes()
		pbundles = partial.get_bundles()
		
		# build partial node lookup.
		pnlookup = dict()
		for i in range(pnodes.size):
			pnlookup[ pnodes[i]['idx'] ] = i
		
		# build set of nodes from partial.
		pnset = set( list( pnodes[:]['idx'] ) )
		
		# flip nodes.
		o_remove = list()
		for i in range(pnodes.size):
			pnodes[i]['orien'] = 1 - pnodes[i]['orien']

		# verify bundles.
		for i in range(pbundles.size):
		
			# simplify.
			X = pbundles[i]['X']
			Y = pbundles[i]['Y']
			
			# flip if any is 1.
			if X == 1 or Y == 1:
				pbundles[i]['X'] = 1 - pbundles[i]['X']
				pbundles[i]['Y'] = 1 - pbundles[i]['Y']
				
		# apply new solution.
		status = self.apply_partial(partial, tri=True)
		if status == False:
			self._just_nodes(pnodes)
			#self._just_bundles(pbundles)
			x = 1
			
			
		# add nodes to tracker.
		self._nodes_added = self._nodes_added.union(pnset)			
	
	
	def apply_partial(self, partial, tri=False):
		''' applies a partial solution to the complete '''
		
		# get solutions.
		pnodes = partial.get_nodes()
		pbundles = partial.get_bundles()
		
		# build partial node lookup.
		pnlookup = dict()
		for i in range(pnodes.size):
			pnlookup[ pnodes[i]['idx'] ] = i
		
		# build set of nodes from partial.
		pnset = set( list( pnodes[:]['idx'] ) )
		
		# check intersection.
		isec = self._nodes_added.intersection(pnset)
		isz = len(isec)
		
		# no-overlap.
		if isz == 0:
			
			# just add nodes.
			self._just_nodes(pnodes)
					
			# just add bundles.
			self._just_bundles(pbundles)
		
		# must merge solution.	
		else:

			# verify nodes.
			for n in isec:
				
				# simplify.
				o_og = self._sol_nodes[n]['orien']
				o_new = pnodes[pnlookup[n]]['orien']
				
				# check orientation the same.
				if o_og != o_new:
					if tri == True:
						return False
					else:
						logging.error("error in bi-comp merging of orientation")
						sys.exit(1)
			
			# verify bundles.
			for i in range(pbundles.size):
			
				# get index.
				idxa = pbundles[i]['idxa']
				idxb = pbundles[i]['idxb']
				idxbun = self._bundle_lookup[self._bundle_key(idxa, idxb)]
				
				# skip new.
				if idxa not in isec or idxb not in isec: continue
				
				# simplify.
				X_og = self._sol_bundles[idxbun]['X']
				Y_og = self._sol_bundles[idxbun]['Y']
				X_new = pbundles[i]['X']
				Y_new = pbundles[i]['Y']
				
				# verify path variables are consistent.
				if X_og != X_new or Y_og != Y_new:
					if tri == True:
						return False
					else:
						print pbundles[i]
						print self._sol_bundles[idxbun]
						logging.error("error in bi-comp merging of paths")
						sys.exit(1)
				
			# now just apply nodes.
			self._just_nodes(pnodes)
				
			# now just apply bundles.
			self._just_bundles(pbundles)
			
		# add nodes to tracker.
		self._nodes_added = self._nodes_added.union(pnset)
		
		# return it worked.
		return True
		
	def apply_agp(self, bundles, agp_file):
		''' applies an AGP solution '''
		
		# sanity check.
		if self._sol_added > 0:
			logging.error("can't apply AGP after solutions added")
			sys.exit(1)
		
		# load the agp array.
		agp_edges = load_agps(agp_file)
		
		# ensure sorted by scaffname and scafidx.
		agp_edges.sort(order=['scaf_name','scaf_idx'])
		
		# apply orientation solutions.
		for i in range(agp_edges.size):
			
			# skip non contigs.
			if agp_edges[i]['comp_type'] != "W": continue
			
			# lookup index.
			idxa = self._nindex(agp_edges[i]['comp_name'])
			
			# apply orientation.
			self._sol_nodes[idxa]['idx'] = idxa
			self._sol_nodes[idxa]['orien'] = agp_edges[i]['comp_orien']
			
			# add to added var.
			self._nodes_added.add(idxa)
			
		# build list of component offsets.
		offsets = dict()
		for i in range(agp_edges.size):

			# skip non contigs.
			if agp_edges[i]['comp_type'] != "W": continue
			
			# record index.
			if agp_edges[i]['scaf_name'] not in offsets:
				offsets[agp_edges[i]['scaf_name']] = list()
			offsets[agp_edges[i]['scaf_name']].append(i)
					
			
		# grow bundle array by this size.
		to_grow = 0
		idxbun = self._sol_bundles.size
		for key in offsets:
			to_grow += len(offsets[key]) - 1
		self._sol_bundles.resize(idxbun + to_grow)
		
		# add bundle info to this.
		gaps = dict()
		for key in offsets:
			
			# loop over edges.
			for i in range(len(offsets[key]) - 1):
				
				# get AGP edge.
				ea = agp_edges[offsets[key][i]]
				eb = agp_edges[offsets[key][i+1]]
				
				# get index.
				idxa = self._nindex(ea['comp_name'])
				idxb = self._nindex(eb['comp_name'])
				
				# get gap.
				gaps[(idxa,idxb)] = eb['scaf_start'] - ea['scaf_stop']
				
				# add to bundles.
				self._sol_bundles[idxbun]['idxa'] = idxa
				self._sol_bundles[idxbun]['idxb'] = idxb
				self._sol_bundles[idxbun]['X'] = 1
				idxbun += 1
				
		# default the state variables.
		self._sol_bundles[:]['S'] = -1
		self._sol_bundles[:]['A'] = -1
		self._sol_bundles[:]['B'] = -1
		self._sol_bundles[:]['C'] = -1
		self._sol_bundles[:]['D'] = -1

		# return the gap estimates.
		return gaps
		

	def _nindex(self, a):
		''' return the index in bundle solution array of the pair '''
		return self._node_lookup[a]
		
	def _bindex(self, a, b):
		''' return the index in bundle solution array of the pair '''
		if a > b:
			return self._bundle_lookup[(a,b)]
		else:
			return self._bundle_lookup[(b,a)]
