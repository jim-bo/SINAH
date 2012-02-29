'''
provides function to estimate gap size by average.
'''

def _est_dist(nodea, nodeb, edge):
	''' estimates the distance '''
	
	# simplify.
	ctg_a_sz = nodea['ctg_width']
	ctg_b_sz = nodeb['ctg_width']
	
	aid = edge['ctg_a_idx']
	bid = edge['ctg_b_idx']
	ao = edge['read_a_orien']
	bo = edge['read_b_orien']
	al = edge['read_a_left_pos']
	ar = edge['read_a_right_pos']
	bl = edge['read_b_left_pos']
	br = edge['read_b_right_pos']
	
	# arrange contig size.
	if aid < bid:
		asz = ctg_a_sz
		bsz = ctg_b_sz
	else:
		asz = ctg_b_sz
		bsz = ctg_a_sz					

	# is this type one.
	da = db = None
	if aid < bid:
		
		# determine happyness class.
		if ao == 0 and bo == 0:
			da = asz - ar
			db = bl
		elif ao == 0 and bo == 1:
			da = asz
			db = bsz - br
		elif ao == 1 and bo == 0:
			da = al
			db = bl
		elif ao == 1 and bo == 1:
			da = al
			db = bsz - br
		
	# no its type two.
	else:
		# determine happyness class.
		if ao == 1 and bo == 1:
			da = al
			db = bsz - br
		elif ao == 0 and bo == 1:
			da = asz - ar
			db = bsz - br
		elif ao == 1 and bo == 0:
			da = al
			db = bl
		elif ao == 0 and bo == 0:
			da = asz - ar
			db = bl

	# return it
	return float(edge['insert_size'] - da - db)
	

def gap_avg(nodes, edges, bundles, agps):
	''' computes gap estimate '''
			
	# setup gap estimates.
	gap_estimate = dict()
	for i in range(bundles.size):
		
		# make keys.
		idxa = bundles[i]['ctg_a_idx']
		idxb = bundles[i]['ctg_b_idx']
		if idxa < idxb:
			key = (idxa, idxb)
		else:
			key = (idxb, idxa)
		
		# add to dictionary.
		gap_estimate[key] = {'count':0.0, 'total':0.0}
			
	# identifies edges for each bundle.
	for i in range(edges.size):
		
		# skip invalid.
		if edges[i]['invalid'] == True:
			continue
		
		# get key.
		idxa = edges[i]['ctg_a_idx']
		idxb = edges[i]['ctg_b_idx']
		if idxa < idxb:
			key = (idxa, idxb)
		else:
			key = (idxb, idxa)		
		
		# sanity check.
		if key not in gap_estimate:
			continue
		
		# call dist estimate function.
		dist = _est_dist(nodes[idxa], nodes[idxb], edges[i])
			
		# track it.
		gap_estimate[key]['total'] += dist
		gap_estimate[key]['count'] += 1.0
		
	# finalize it.
	for key in gap_estimate:
					
		# normal.
		if gap_estimate[key]['count'] != 0.0:
			gap_estimate[key] = gap_estimate[key]['total'] / gap_estimate[key]['count']
		else:
			
			# look in agps.
			if key in agps:
				gap_estimate[key] = agps[key]
			elif (key[1], key[0]) in agps:
				gap_estimate[key] = agps[(key[1], key[0])]
			else:
				logging.error("gap errrrrrror")
				sys.exit(1)
			
	
	# return the estimates.
	return gap_estimate
