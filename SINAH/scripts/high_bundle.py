'''
plot true positive bundle size vs false positive bundle size
'''

# program imports.
from data_structs.nodes import load_nodes
from data_structs.nodes import create_lookup
from data_structs.edges import load_edges
from data_structs.bundles import load_bundles
from data_structs.agps import load_agps
from data_structs import types

# system imports.
import time
import random
import itertools
import networkx as nx
import logging
import sys
import os
from scitools.std import *

time_start = time.time()

# logging.
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )
#logging.basicConfig(level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s', )

# parameters.
input_nodes_file = os.path.abspath(sys.argv[1])
input_edges_file = os.path.abspath(sys.argv[2])
input_agp_file = os.path.abspath(sys.argv[3])

#################### functions ####################

def draw(name, subg):
    ''' draw an interactive graph.'''
    # write dot file.
    dot = "./tmp.dot"
    pic = "/home/jlindsay/central/transfer/%s.jpg" % name
    
    nx.write_dot(subg, dot)
    
    # convert to picture using neato.
    subprocess.call(["neato", "-Tjpg", "-o", pic, dot])
    
    # remove dot.
    subprocess.call(["rm","-f",dot])

def create_dir(file_path):
    ''' creates a dir. '''
    if os.path.isdir(file_path) == False:
        subprocess.call(["mkdir",file_path])

def make_key(a, b):
    ''' makes sorted key'''
    if a < b:
        return (a,b)
    else:
        return (b,a)

########### script ################## 

# load hdf5 information.
logging.info("loading data arrays")
nodes = load_nodes(input_nodes_file)
edges = load_edges(input_edges_file)
agps = load_agps(input_agp_file)

nlookup = create_lookup(nodes)

# build bundle count.
logging.info("counting bundles")
blookup = dict()
for i in range(edges.size):

    # get id
    idxa = edges[i]['ctg_a_idx']
    idxb = edges[i]['ctg_b_idx']
    key = make_key(idxa, idxb)
    
    # count it.
    if key not in blookup:
        blookup[key] = 0
    blookup[key] += 1
    
# build adjacency set.
logging.info("building true positive set")
tp = set()
for i in range(agps.size):
    
    # skip to gaps.
    if agps[i]['comp_type'] == 'N':

        # get id
        idxa = agps[i-1]['comp_name']
        idxb = agps[i+1]['comp_name']
        idxa = nlookup[idxa]
        idxb = nlookup[idxb]
        key = make_key(idxa, idxb)        
        
        # add to set.
        tp.add(key)
        
# evaluate bundles.
glist = []
blist = []
largest = 0
for key in blookup:
    
    # get largest.
    if blookup[key] > largest:
        largest = blookup[key]
    
    # check if true or not.
    if key in tp:
        glist.append(blookup[key])
    else:
        blist.append(blookup[key])
        
# computer values.
x = [i for i in range(0,largest+1)]
y1 = [0 for i in range(0,largest+1)]
y2 = [0 for i in range(0,largest+1)]
for a in glist:
    y1[a] += 1
for a in blist:
    y2[a] += 1
        
# plot these.
plot(x, y1, 'b', x, y2, 'r', xlabel='bundle size', ylabel='count',
     legend=('blue=TP', 'red=FP'),
     title='Assesing Bundle Size as a Filtering Parameter: 1%',
     hardcopy='bsize_filter.eps'
)
