'''
Created on Mar 21, 2011

solves ILP without path constraints

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
	('WT_A', np.double),\
	('WT_B', np.double),\
	('WT_C', np.double),\
	('WT_D', np.double)\
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
		self._neibs = dict()
		
		# add the constraints.
		self._constrain()
		
		# set loaded.
		self._loaded = True		
			
		
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
		del self._neibs
		self._cpx = None
		
	
		# clear loaded.
		self._loaded = False
		self._solved = False
		
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
		
	def force_one(self, n, val):
		''' force orientation of 1 node '''
		
		# sanity check 1.
		if self._loaded != True:
			logging.error("can't add bi-comp constraints to non loaded ILP")
			sys.exit(1)
		
		# build vars.
		Si = "S_%i" % n
			
		# sanity check 2.
		if Si not in self._var_defined:
			logging.error("error in bicomp constraint: %s" % Si)
			sys.exit(1)
			
		# make constraint.
		c1 = cplex.SparsePair(ind=[Si], val=[1])
		self._cpx.linear_constraints.add( lin_expr = [c1], senses = ["E"], rhs = [val], names = ["force1"] )
		
	def force_path(self, a, b):
		''' force orientation of 1 node '''
		
		# sanity check 1.
		if self._loaded != True:
			logging.error("can't add bi-comp constraints to non loaded ILP")
			sys.exit(1)
		
		# build vars.
		Xij = "X_%i_%i" % (a, b)
		Xji = "X_%i_%i" % (b, a)
			
		# sanity check 2.
		if Xij not in self._var_defined or Xji not in self._var_defined:
			logging.error("error in bicomp constraint: %s %s" % (Xij,Xji))
			sys.exit(1)
			
		# make constraint.
		c1 = cplex.SparsePair(ind=[Xij], val=[1])
		c2 = cplex.SparsePair(ind=[Xji], val=[1])
		self._cpx.linear_constraints.add( lin_expr = [c1, c2], senses = ["E", "E"], rhs = [1, 0], names = ["force2", "force2"] )

	def weight_one(self, n, val):
		''' weight orientation of 1 node '''
		
		# sanity check 1.
		if self._loaded != True:
			logging.error("can't add bi-comp constraints to non loaded ILP")
			sys.exit(1)
		
		# build vars.
		Si = "S_%i" % n
			
		# sanity check 2.
		if Si not in self._var_defined:
			logging.error("error in bicomp constraint: %s" % Si)
			sys.exit(1)
			
		# make constraint.
		self._cpx.objective.set_linear(Si, val)

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
			
			# get objective value.
			obj_val = self._cpx.solution.get_objective_value()
			
		except CplexSolverError, e:
			
			# if no solution found return empty sol and -1.
			self._sol = None
			obj_val = -1
			
			#logging.error("exception raised during solve: " + str(e))
			#sys.exit(1)
	
		# set solved to true.
		self._solved = True
		
		# return solution.
		return self._sol, obj_val
		
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
			
			'''
			if val_Aij + val_Bij + val_Cij + val_Dij != 1:
				logging.error("error in orientation solution: %i %i" % (idxa, idxb))
				print val_Sij, val_Aij, val_Bij, val_Cij, val_Dij, val_Xij, val_Xji
				print self._blist[i]
				sys.exit(1)
			'''
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
		
		
	def _constrain_xvars(self, supplied=False):
		''' ensure we have a path using Xij,Xji'''

		# switch on source.
		if supplied == False:
			size = self._blist.size
			active = self._blist
		else:
			size = len(supplied)
			active = supplied

		# loop over bundles.
		for i in range(size):
			
			# simplify.
			idxa = active[i]['idxa']
			idxb = active[i]['idxb']
			
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
		
	def _constrain_two(self, supplied=False):
		''' ensures no two cycles '''

		# switch on source.
		if supplied == False:
			size = self._blist.size
			active = self._blist
		else:
			size = len(supplied)
			active = supplied
		

		# loop over bundles.
		for i in range(size):
			
			# simplify.
			idxa = active[i]['idxa']
			idxb = active[i]['idxb']
			
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
			c4 = cplex.SparsePair(ind=[Aij, Bij, Cij, Dij     ], val=[ 1, 1, 1, 1   ])
		
			# add to ILP.
			self._cpx.linear_constraints.add(\
				lin_expr = [c1, c2, c3, c4],\
				senses = ["L", "L", "L", "E"],\
				rhs = [1, 1, 0, 1],\
				names = ["cycle2", "cycle2", "cycle2", "cycle2"]\
			)

		
	def _constrain_sij(self, supplied=False):
		''' fixes Sij = 0 when Si=Sj, 1 otherwise.'''
		
		# switch on source.
		if supplied == False:
			size = self._blist.size
			active = self._blist
		else:
			size = len(supplied)
			active = supplied
		
		# loop over bundles.
		for i in range(size):
			
			# simplify.
			idxa = active[i]['idxa']
			idxb = active[i]['idxb']
			
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
			
			# build cooefficients.
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

class BiConstraints(object):
	''' encapsulates bi-constraints used in bottom up '''
	
	def __init__(self):
		''' init '''
		
		# perpare path var lists.
		self._pvars = list()
		
		
	def load(self, p, bisols):
		''' adds these 4 constraints '''
		
		# create path variables.
		Xar = "X_%i_r" % p
		Xra = "X_r_%i" % p
		
		Xat = "X_%i_t" % p
		Xta = "X_t_%i" % p
		
		# add to list.
		self._pvars.append(Xar)
		self._pvars.append(Xra)
		self._pvars.append(Xat)
		self._pvars.append(Xta)
		
		
	
