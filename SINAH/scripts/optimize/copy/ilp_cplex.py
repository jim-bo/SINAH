'''
Created on Mar 21, 2011

@author: jlindsay
'''

# system imports.
import sys
import logging
import numpy as np
import cplex
import math
from cplex.exceptions import CplexSolverError

nlist_dt = np.dtype([\
	('idx', np.int)\
])

blist_dt = np.dtype([\
	('idxa', np.int),\
	('idxb', np.int),\
	('WT_A', np.int),\
	('WT_B', np.int),\
	('WT_C', np.int),\
	('WT_D', np.int)\
])

tlist_dt = np.dtype([\
	('idxa', np.int),\
	('idxb', np.int),\
	('idxc', np.int)\
])
	

class SpqrIlp(object):
	'''
	implements SPQR ILP using CPLEX
	'''

	def __init__(self, log_file, err_file):
		'''
		constructor
		'''	
		
		# save file ptrs.
		self._log_file = log_file
		self._err_file = err_file
		
		# clear logs.
		tmp = open(self._log_file, "w")
		tmp.close()
		tmp = open(self._err_file, "w")
		tmp.close()
		
		# set loaded var.
		self._loaded = False
		self._solved = False
		

	def load(self, nlist, blist, tlist, sol):
		''' loads the node, bundle lists then creates triangle list'''
		
		# sanity check.
		if self._loaded == True:
			logging.error("ILP already loaded.")
			sys.exit(1)
		
		# save pointers.
		self._nlist = nlist
		self._blist = blist
		self._tlist = tlist
		self._sol = sol
		
		# initiate cplex object.
		self._cpx = cplex.Cplex()
		
		# set log files.
		self._cpx.set_log_stream(self._log_file)
		self._cpx.set_results_stream(self._log_file)
		self._cpx.set_warning_stream(self._err_file)
		self._cpx.set_error_stream(self._err_file)
		
		# prepare lookup structures.
		self._var_defined = set()
		self._node_lookup = dict()
		self._pair_lookup = dict()
		
		# add the constraints.
		self._constrain()
		
		# set loaded.
		self._loaded = True		
			
	def triclear(self):
		''' remove tri-forced constraint '''
		
		# remove tri-forced constraint.
		self._cpx.linear_constraints.delete("triforce")		
	
		# clear the solved key.
		self._solved = False
	
	def changesol(self, sol):
		''' changes solution object '''
		self._sol = sol
		
	def clear(self):
		''' resets ILP completely '''
		
		# sanity.
		if self._cpx == None:
			logging.error("ILP already deleted")
			sys.exit(1)
			
		# sanity.
		if self._solved == False:
			logging.error("ILP not solved")
			sys.exit(1)
		
		# remove cplex and other vars.
		del self._cpx
		del self._nlist
		del self._blist
		del self._tlist
		del self._var_defined
		self._cpx = None
		
	
		# clear loaded.
		self._loaded = False
		self._solved = False

	def trimod(self, constraints):
		''' adds tri-mod constraints '''
		
		# sanity check 1.
		if self._loaded != True:
			logging.error("can't add tri-mod constraints to non loaded ILP")
			sys.exit(1)
			
		# loop over contraints.
		for cut, vals in constraints:
			
			# prep var.
			Sij = "S_%i_%i" % (cut[0], cut[1])
			
			# add if not there.
			if Sij not in self._var_defined:
				self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Sij] )
			
			
			# simplify.
			val_s = vals[0]
			val_d = vals[1]
			
			# add constraints.
			self._mod_objective(Sij, cut[0], cut[1], val_s, val_d)
			
			
	
	def tricomp(self, cut_x, cut_y, val, vals):
		''' adds tri-comp constraints '''
		
		# sanity check 1.
		if self._loaded != True:
			logging.error("can't add tri-comp constraints to non loaded ILP")
			sys.exit(1)
			
		# make vars.
		Sx = "S_%i" % cut_x
		Sy = "S_%i" % cut_y

		# sanity check 2.
		if Sx not in self._var_defined and Sy not in self._var_defined:
			logging.error("tri-constraints error 2: %s %s" % (Sx, Sy))
			sys.exit(1)
			
		# add constraint.
		self._constrain_tri(cut_x, cut_y, val, vals)
		
	def _constrain_tri(self, cut_x, cut_y, val, vals):
		''' adds tri-comp same/diff constraints '''
		
		# make vars.
		Sx = "S_%i" % cut_x
		Sy = "S_%i" % cut_y
		
		# build constraint.
		if val == 0:
			c1 = cplex.SparsePair( ind = [Sx, Sy], val = [1, -1] )
			self._cpx.linear_constraints.add( \
				lin_expr = [c1],\
				senses = ["E"],\
				rhs = [0],\
				names = ["triforce"] 
			)
		else:
			c1 = cplex.SparsePair( ind = [Sx, Sy], val = [1, 1] )
			self._cpx.linear_constraints.add( \
				lin_expr = [c1],\
				senses = ["G"],\
				rhs = [0],\
				names = ["triforce"] 
			)
		
		# constrain values.
		if vals != False:
			c1 = cplex.SparsePair( ind = [Sx], val = [1] )
			c2 = cplex.SparsePair( ind = [Sy], val = [1] )
			self._cpx.linear_constraints.add( \
				lin_expr = [c1, c2],\
				senses = ["E", "E"],\
				rhs = [vals[0], vals[1]],\
				names = ["triforce", "triforce"] 
			)			
		
	def bicomp(self, active, sol_nodes, sol_pairs, tri=False):
		''' adds bi-comp intersection contraints '''
		
		# sanity check 1.
		if self._loaded != True:
			logging.error("can't add bi-comp constraints to non loaded ILP")
			sys.exit(1)
			
		# sanity check 2.
		for e in sol_nodes:
			
			# build vars.
			Si = "S_%i" % e['idx']
			
			# check we have orientation vars.
			if Si not in self._var_defined:
				if tri == False:
					logging.error("error in bicomp constraint: %s" % Si)
					sys.exit(1)
				else:
					self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Si] )
					self._var_defined.add(Si)
				
		# add vars if necessary.
		for e in sol_pairs:
			
			# simplify.
			idxa = e['idxa']
			idxb = e['idxb']
			
			# build vars.
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)
			
			# check if we need to add a var.
			if Xij not in self._var_defined:
				self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Xij] )
				self._var_defined.add(Xij)

			if Xji not in self._var_defined:
				self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Xji] )
				self._var_defined.add(Xji)
			
		# add constraints.
		self._constrain_bi(sol_nodes, sol_pairs)
	
	def pizza(self, active, sol_pairs):
		''' adds tri-comp intersection contraints '''
		
		# sanity check 1.
		if self._loaded != True:
			logging.error("can't add bi-comp constraints to non loaded ILP")
			sys.exit(1)
		
				
		# add vars if necessary.
		for e in sol_pairs:
			
			# simplify.
			idxa = e['idxa']
			idxb = e['idxb']
			
			# build vars.
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)
			
			# check if we need to add a var.
			if Xij not in self._var_defined:
				self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Xij] )
				self._var_defined.add(Xij)

			if Xji not in self._var_defined:
				self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Xji] )
				self._var_defined.add(Xji)
			
		# add constraints.
		self._constrain_pizza(sol_pairs)
	


	
	def solve(self, file_path=False):
		''' runs the ilp on loaded info '''

		# sanity check.
		if self._loaded == False:
			logging.error("ILP not loaded.")
			sys.exit(1)			
			
		# sanity check.
		if self._solved == True:
			logging.error("shouldn't solve ILP twice.")
			sys.exit(1)			
			
		# write ILP to file.
		if file_path != False:
			self._cpx.write(file_path, filetype="lp")
			
		# call the solve code.
		try:
			
			# call the solve method.	
			self._cpx.solve()
			
			# populate solution.
			self._populate_sol()
			
		except CplexSolverError, e:
			logging.error("exception raised during solve: " + str(e))
			sys.exit(1)
	
		# set solved to true.
		self._solved = True
		
		# return solution.
		return self._sol, self._cpx.solution.get_objective_value()
		
	def _populate_sol(self):
		''' populates solution object after running '''
		
		# loop over nodes.
		for i in range(self._nlist.size):
			
			# build var.
			idx = self._nlist[i]['idx']
			Si = "S_%i" % idx
			
			# get result.
			val = int(self._cpx.solution.get_values(Si))
			
			# add to solution.
			self._sol.set_node(idx, val)
			
		# loop over bundles.
		for i in range(self._blist.size):
			
			# simplify.
			idxa = self._blist[i]['idxa']
			idxb = self._blist[i]['idxb']
			
			# build vars.
			Sij = "S_%i_%i" % (idxa, idxb)
			Aij = "A_%i_%i" % (idxa, idxb)
			Bij = "B_%i_%i" % (idxa, idxb)
			Cij = "C_%i_%i" % (idxa, idxb)
			Dij = "D_%i_%i" % (idxa, idxb)
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)		
			
			# get values.
			val_Sij = int(self._cpx.solution.get_values(Si))
			val_Aij = int(self._cpx.solution.get_values(Aij))
			val_Bij = int(self._cpx.solution.get_values(Bij))
			val_Cij = int(self._cpx.solution.get_values(Cij))
			val_Dij = int(self._cpx.solution.get_values(Dij))
			val_Xij = int(self._cpx.solution.get_values(Xij))
			val_Xji = int(self._cpx.solution.get_values(Xji))
			
			# sanity check.
			if val_Xij == val_Xji and val_Xij == 1:
				logging.error("error in path logic: %s %s" % (Xij, Xji))
				sys.exit(1)
			
			# add to solution.
			self._sol.set_bundle(idxa, idxb, val_Sij, val_Aij, val_Bij, val_Cij, val_Dij, val_Xij, val_Xji)
		
		
	
		
	def _constrain(self):
		''' adds constraints to CPLEX '''
		
		# add node variables.
		self._add_node_vars()
		
		# add pair variables.
		self._add_pair_vars()
		
		# set the objective.
		self._set_objective()
		
		# add the constraints to fix Sij behavior.
		self._constrain_sij()
		
		# ensure no 2 cycles.
		self._constrain_two()
		
		# ensure no 3 cycles.
		self._constrain_three()
		
		# ensure Xij, Xji function properly.
		self._constrain_xvars()
		
		# ensure we calculate a path.
		self._constrain_path()
		
		
	def _constrain_path(self):
		''' ensures ILP calculates a path '''
		
		# build list of neighbors.
		self._neibs = dict()
		for i in range(self._blist.size):
			
			# bootstrap.
			if self._blist[i]['idxa'] not in self._neibs:
				self._neibs[self._blist[i]['idxa']] = set()
			if self._blist[i]['idxb'] not in self._neibs:
				self._neibs[self._blist[i]['idxb']] = set()

			# track.
			self._neibs[self._blist[i]['idxa']].add(self._blist[i]['idxb'])
			self._neibs[self._blist[i]['idxb']].add(self._blist[i]['idxa'])
			
		# loop over nodes.
		for idxa in self._nlist[:]['idx']:
			
			# neighbor check.
			if idxa not in self._neibs:
				continue
			
			# loop over each neighbor.
			in_var = []
			in_val = []
			out_var = []
			out_val = []
			for idxb in self._neibs[idxa]:
				
				# set vars.
				Xij = "X_%i_%i" % (idxa, idxb)
				Xji = "X_%i_%i" % (idxb, idxa)
				
				# skip if not defined.
				if Xij not in self._var_defined or Xji not in self._var_defined:
					continue
				
				# outgoing.
				out_var.append(Xij)
				out_val.append(1)
				
				# incoming.
				in_var.append(Xji)
				in_val.append(1)
						
			# create and add constraints.
			c1 = cplex.SparsePair( ind = in_var, val = in_val )
			c2 = cplex.SparsePair( ind = out_var, val = out_val )
			self._cpx.linear_constraints.add(
				lin_expr = [c1, c2],\
				senses = ["L", "L"],\
				rhs = [1, 1],\
				names = ["pathsum", "pathsum"]\
			)

	def _constrain_pizza(self, sol_pairs):
		''' adds constraints to force Si, Xij, Xji be same as already solved for '''
		'''
		# copy subset of neighbor dictionary.
		neibs2 = dict()
		for e in sol_nodes:
			if e['idx'] in self._neibs:
				neibs2[e['idx']] = self._neibs[e['idx']]
		
		# add AGP inspired neighbors.
		for e in sol_pairs:
			
			# simplify.
			idxa = e['idxa']
			idxb = e['idxb']
			
			# build neighbors.
			if idxa not in neibs2:
				neibs2[idxa] = set()
			neibs2[idxa].add(idxb)
				
			if idxb not in neibs2:
				neibs2[idxb] = set()
			neibs2[idxb].add(idxa)		
		'''
		# force path.
		for e in sol_pairs:
			
			# simplify.
			idxa = e['idxa']
			idxb = e['idxb']
			pathij = e['X']
			pathji = e['Y']
			
			# build vars.
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)
			
			print "hamster", Xij, Xji, pathij, pathji
			
			# build constraint.
			c1 = cplex.SparsePair( ind = [Xij], val = [1] )
			c2 = cplex.SparsePair( ind = [Xji], val = [1] )
			self._cpx.linear_constraints.add( \
				lin_expr = [c1, c2],\
				senses = ["E", "E"],\
				rhs = [pathij, pathji],\
				names = ["force_path", "force_path"] 
			)
			'''
			# update new neighbor list.
			if idxa in neibs2:
				neibs2[idxa].add(idxb)
			if idxb in neibs2:
				neibs2[idxb].add(idxa)
			'''
		'''
		# update path sum variables.
		for idxa in neibs2:
			
			# loop over each neighbor.
			in_var = []
			in_val = []
			out_var = []
			out_val = []
			for idxb in neibs2[idxa]:
				
				# set vars.
				Xij = "X_%i_%i" % (idxa, idxb)
				Xji = "X_%i_%i" % (idxb, idxa)
				
				# skip if not defined.
				if Xij not in self._var_defined or Xji not in self._var_defined:
					logging.error("error in bi-comp path update")
					sys.exit(1)
				
				# outgoing.
				out_var.append(Xij)
				out_val.append(1)
				
				# incoming.
				in_var.append(Xji)
				in_val.append(1)
						
			# create and add constraints.
			c1 = cplex.SparsePair( ind = in_var, val = in_val )
			c2 = cplex.SparsePair( ind = out_var, val = out_val )
			self._cpx.linear_constraints.add(
				lin_expr = [c1, c2],\
				senses = ["L", "L"],\
				rhs = [1, 1],\
				names = ["pathsum_update", "pathsum_update"]\
			)		
		'''
	def _constrain_bi(self, sol_nodes, sol_pairs):
		''' adds constraints to force Si, Xij, Xji be same as already solved for '''
		
		# force orientation.
		for e in sol_nodes:
			
			# simplify.
			idx = e['idx']
			orien = e['orien']
			
			# build vars.
			Si = "S_%i" % e['idx']
			
			# build constraint.
			c1 = cplex.SparsePair( ind = [Si], val = [1] )
			self._cpx.linear_constraints.add( \
				lin_expr = [c1],\
				senses = ["E"],\
				rhs = [orien],\
				names = ["force_orien"] 
			)				

		# copy subset of neighbor dictionary.
		neibs2 = dict()
		for e in sol_nodes:
			if e['idx'] in self._neibs:
				neibs2[e['idx']] = self._neibs[e['idx']]
		
		# add AGP inspired neighbors.
		for e in sol_pairs:
			
			# simplify.
			idxa = e['idxa']
			idxb = e['idxb']
			
			# build neighbors.
			if idxa not in neibs2:
				neibs2[idxa] = set()
			neibs2[idxa].add(idxb)
				
			if idxb not in neibs2:
				neibs2[idxb] = set()
			neibs2[idxb].add(idxa)		
		
		# force path.
		for e in sol_pairs:
			
			# simplify.
			idxa = e['idxa']
			idxb = e['idxb']
			pathij = e['X']
			pathji = e['Y']
			
			# build vars.
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)
			
			# build constraint.
			c1 = cplex.SparsePair( ind = [Xij], val = [1] )
			c2 = cplex.SparsePair( ind = [Xji], val = [1] )
			self._cpx.linear_constraints.add( \
				lin_expr = [c1, c2],\
				senses = ["E", "E"],\
				rhs = [pathij, pathji],\
				names = ["force_path", "force_path"] 
			)

			# update new neighbor list.
			if idxa in neibs2:
				neibs2[idxa].add(idxb)
			if idxb in neibs2:
				neibs2[idxb].add(idxa)

			
		# update path sum variables.
		for idxa in neibs2:
			
			# loop over each neighbor.
			in_var = []
			in_val = []
			out_var = []
			out_val = []
			for idxb in neibs2[idxa]:
				
				# set vars.
				Xij = "X_%i_%i" % (idxa, idxb)
				Xji = "X_%i_%i" % (idxb, idxa)
				
				# skip if not defined.
				if Xij not in self._var_defined or Xji not in self._var_defined:
					logging.error("error in bi-comp path update")
					sys.exit(1)
				
				# outgoing.
				out_var.append(Xij)
				out_val.append(1)
				
				# incoming.
				in_var.append(Xji)
				in_val.append(1)
						
			# create and add constraints.
			c1 = cplex.SparsePair( ind = in_var, val = in_val )
			c2 = cplex.SparsePair( ind = out_var, val = out_val )
			self._cpx.linear_constraints.add(
				lin_expr = [c1, c2],\
				senses = ["L", "L"],\
				rhs = [1, 1],\
				names = ["pathsum_update", "pathsum_update"]\
			)		
		
	def _constrain_three(self):
		''' ensure no 3 cycles '''
	
		# loop over triangles.
		for i in range(self._tlist.size):
			
			# simplify.
			idxa = self._tlist[i]['idxa']
			idxb = self._tlist[i]['idxb']
			idxc = self._tlist[i]['idxc']
	
			# make variables.
			Xij = "X_%i_%i" % (idxa, idxb)
			Xjk = "X_%i_%i" % (idxb, idxc)
			Xki = "X_%i_%i" % (idxc, idxa)

			Xji = "X_%i_%i" % (idxb, idxa)
			Xkj = "X_%i_%i" % (idxc, idxb)
			Xik = "X_%i_%i" % (idxa, idxc)
			
			# sanity check.
			for x in [Xij, Xjk, Xki, Xji, Xkj, Xik]:
				
				# check if variable is there.
				if x not in self._var_defined:
					logging.error("bad logic in triangles, fix this programmer!: %s" % x)	
					sys.exit(1)
			
			# add two triangle constraints.
			c1 = cplex.SparsePair( ind = [Xij, Xjk, Xki], val = [1,1,1] )
			c2 = cplex.SparsePair( ind = [Xji, Xkj, Xik], val = [1,1,1] )
			self._cpx.linear_constraints.add( lin_expr = [c1, c2], senses = ["L", "L"], rhs = [2, 2] )	
		
		
	def _constrain_xvars(self):
		''' ensure we have a path using Xij,Xji'''

		# loop over bundles.
		for i in range(self._blist.size):
			
			# simplify.
			idxa = self._blist[i]['idxa']
			idxb = self._blist[i]['idxb']
			
			# build vars.
			Si = "S_%i" % idxa
			Sj = "S_%i" % idxb
			Sij = "S_%i_%i" % (idxa, idxb)
			Aij = "A_%i_%i" % (idxa, idxb)
			Bij = "B_%i_%i" % (idxa, idxb)
			Cij = "C_%i_%i" % (idxa, idxb)
			Dij = "D_%i_%i" % (idxa, idxb)
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)

			# create and add Xij constraints.
			c1 = cplex.SparsePair( ind = [Aij, Si, Sj, Xij], val = [2,-1,-1,-2] )
			c2 = cplex.SparsePair( ind = [Dij, Si, Sj, Xij], val = [2, 1, 1,-2] )
			c3 = cplex.SparsePair( ind = [Bij,     Sj, Xij], val = [1,    1,-1] )
			c4 = cplex.SparsePair( ind = [Cij, Si,     Xij], val = [1, 1,   -1] )
			names = ["xvars_ij", "xvars_ij", "xvars_ij", "xvars_ij"]

			self._cpx.linear_constraints.add( \
				lin_expr = [c1, c2, c3, c4],\
				senses = ["L", "L", "L", "L"],\
				rhs = [0, 2, 1, 1],\
				names = names\
			)		
			
			# create and add Xji constraints.
			c1 = cplex.SparsePair( ind = [Aij, Si, Sj, Xji], val = [2, 1, 1,-2] )
			c2 = cplex.SparsePair( ind = [Dij, Si, Sj, Xji], val = [2,-1,-1,-2] )
			c3 = cplex.SparsePair( ind = [Bij, Si,     Xji], val = [1, 1,   -1] )
			c4 = cplex.SparsePair( ind = [Cij,     Sj, Xji], val = [1,    1,-1] )
			names = ["xvars_ji", "xvars_ji", "xvars_ji", "xvars_ji"]

			self._cpx.linear_constraints.add( \
				lin_expr = [c1, c2, c3, c4],\
				senses = ["L", "L", "L", "L"],\
				rhs = [2, 0, 1, 1],\
				names = names\
			)
		
	def _constrain_two(self):
		''' ensures no two cycles '''

		# loop over bundles.
		for i in range(self._blist.size):
			
			# simplify.
			idxa = self._blist[i]['idxa']
			idxb = self._blist[i]['idxb']
			
			# build vars.
			Sij = "S_%i_%i" % (idxa, idxb)
			Aij = "A_%i_%i" % (idxa, idxb)
			Bij = "B_%i_%i" % (idxa, idxb)
			Cij = "C_%i_%i" % (idxa, idxb)
			Dij = "D_%i_%i" % (idxa, idxb)
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)
		
			# make constraint.
			c1 = cplex.SparsePair(ind=[Xij, Xji], val=[ 1, 1])
			c2 = cplex.SparsePair(ind=[Aij,           Dij, Sij], val=[ 1,       1, 1])
			c3 = cplex.SparsePair(ind=[     Bij, Cij,      Sij], val=[    1, 1,   -1])
		
			# add to ILP.
			self._cpx.linear_constraints.add(\
				lin_expr = [c1, c2, c3],\
				senses = ["L", "L", "L"],\
				rhs = [1, 1, 0],\
				names = ["2cycle", "2cycle", "2cycle"]\
			)
		
	def _constrain_sij(self):
		''' fixes Sij = 0 when Si=Sj, 1 otherwise.'''
		
		# loop over bundles.
		for i in range(self._blist.size):
			
			# simplify.
			idxa = self._blist[i]['idxa']
			idxb = self._blist[i]['idxb']
			
			# build vars.
			Si = "S_%i" % idxa
			Sj = "S_%i" % idxb
			Sij = "S_%i_%i" % (idxa, idxb)
	
			# build constraint variables.
			c1 = cplex.SparsePair(ind=[Sij,Si,Sj], val=[ 1,-1,-1])
			c2 = cplex.SparsePair(ind=[Sij,Si,Sj], val=[ 1, 1, 1])
			c3 = cplex.SparsePair(ind=[Sij,Si,Sj], val=[ 1, 1,-1])
			c4 = cplex.SparsePair(ind=[Sij,Si,Sj], val=[ 1,-1, 1])
			names = ["jin4", "jin4", "jin4", "jin4"]
		
			# add to cplex.
			self._cpx.linear_constraints.add(\
				lin_expr = [c1, c2, c3, c4],\
				senses = ["L", "L", "G", "G"],\
				rhs = [0, 2, 0, 0],\
				names = names\
			)
		
	def _mod_objective(self, Sij, i, j, val_s, val_d):
		''' modifies objective for SPQR '''
		
		# modified existing Sij for diff.
		self._cpx.objective.set_linear(Sij, val_d)
		
		# add new variable for same.
		Vij = "V_%i_%i" % (i, j)
		self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Vij] )
		
		# constrain to be opposite of Sij.
		c1 = cplex.SparsePair(ind=[Sij, Vij], val=[ 1,1])
		self._cpx.linear_constraints.add(\
			lin_expr = [c1],\
			senses = ["E"],\
			rhs = [1],\
			names = ['Vforce']\
		)
		
		# set objective.
		self._cpx.objective.set_linear(Sij, val_s)
		
		
	def _set_objective(self):
		''' sets objective '''
		
		# loop over bundles.
		for i in range(self._blist.size):
			
			# simplify.
			idxa = self._blist[i]['idxa']
			idxb = self._blist[i]['idxb']
			
			# build vars.
			Aij = "A_%i_%i" % (idxa, idxb)
			Bij = "B_%i_%i" % (idxa, idxb)
			Cij = "C_%i_%i" % (idxa, idxb)
			Dij = "D_%i_%i" % (idxa, idxb)
			
			# build constants.
			HA = self._blist[i]['WT_A']
			HB = self._blist[i]['WT_B']
			HC = self._blist[i]['WT_C']
			HD = self._blist[i]['WT_D']
			
			# make list.
			tmp1 = [Aij, Bij, Cij, Dij]
			tmp2 = [HA, HB, HC, HD]
			
			# loop over lists.
			for var, value in zip(tmp1, tmp2):
				
				# add to objective.
				self._cpx.objective.set_linear(var, value)
				
		# set objective type.
		self._cpx.objective.set_sense(self._cpx.objective.sense.maximize)
		
	def _add_pair_vars(self):
		''' adds pair variables '''
		
		# loop over bundles.
		for i in range(self._blist.size):
			
			# simplify.
			idxa = self._blist[i]['idxa']
			idxb = self._blist[i]['idxb']
			
			# build vars.
			Sij = "S_%i_%i" % (idxa, idxb)
			Aij = "A_%i_%i" % (idxa, idxb)
			Bij = "B_%i_%i" % (idxa, idxb)
			Cij = "C_%i_%i" % (idxa, idxb)
			Dij = "D_%i_%i" % (idxa, idxb)
			Xij = "X_%i_%i" % (idxa, idxb)
			Xji = "X_%i_%i" % (idxb, idxa)
			
			# make a list.
			tmp = (Sij, Aij, Bij, Cij, Dij, Xij, Xji)
			
			# loop over list.
			for x in tmp:
				
				# add the variables.
				self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [x] )
				
				# populate lookup.
				self._var_defined.add(x)
		
	def _add_node_vars(self):
		''' adds orientation variables '''
		
		# loop over nodes.
		for i in range(self._nlist.size):
			
			# build var.
			Si = "S_%i" % self._nlist[i]['idx']
			
			# add the variable.
			self._cpx.variables.add( lb = [0], ub = [1], types = ["B"], names = [Si] )
			
			# populate lookup.
			self._var_defined.add(Si)
