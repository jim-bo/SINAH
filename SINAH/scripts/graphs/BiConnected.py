"""
  Constructs the block tree of a graph using dfs in O(n+e)

  Implemented by cloning the igraph implementation
"""
__authors__ = "\n".join(['Jesus Cerquides <cerquide@iiia.csic.es>'])

#    Copyright (C) 2004-2010 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

__all__ = [ 'block_tree',
            'block_tree_node_ypages',
            'block_tree_find',
            'bicomponents',
            'is_biconnected'
            ]

import networkx as nx

def block_treeCC(G,bt,components=None,articulations =None):
    """
    Returns the block tree for the undirected graph G. Can use the set
    of bicomponents and the set of articulations received as parameters.
    Otherwise, it computes the set of bicomponents and articulation points
    The returned block tree is a undirected graph that contains two types of
    nodes:
      articulation points: every node of G that is an articulation point is
        a node of its block tree.
      bicomponents: Every biconnected component of G is a node of its block
        tree. This nodes are graphs themselves, containing all the edges in
        the biconnected component.
    An articulation point is linked to a bicomponent if the node appears in
    the bicomponent.
    The function returns also a dictionary that maps each node of G to itself
    (if it is an articulation point) or to its biconnected component otherwise
    """
    if G.number_of_nodes() == 1:
        bt.add_node(G.nodes()[0])
        return
    if components == None or articulations == None:
        components,articulations= bicomponents(G)
    #bt = nx.Graph()
    for art in articulations:
        bt.add_node(art)
    graph_components = {}
    for component in components:
        graph_components[component] = nx.Graph()
    for (u,v) in G.edges_iter():
        for component in components:
            if u in component and v in component:
                graph_components[component].add_edge(u,v)
    for component in components:
        bt.add_node(graph_components[component])
        for art in component & articulations:
            bt.add_edge(graph_components[component],art)

    #return bt

def block_tree(G,components=None,articulations =None):
    #print G.nodes()
    #print G.edges()
    bt = nx.Graph()
    for CC in nx.algorithms.components.connected.connected_component_subgraphs(G):
        block_treeCC(CC,bt,components,articulations)
    return bt

def block_tree_node_ypages(bt):
    node_ypages = {}
    for node in bt:
        if type(node)==nx.Graph:
            for u in node:
                node_ypages[u] = node
    for node in bt:
        if type(node)!=nx.Graph:
            node_ypages[node] = node
    return node_ypages

##def block_tree_articulationpairs_ypages(bt):
##    for component in components:
##        bt.add_node(graph_components[component])
##        for node in component:
##            if node in articulations:
##                bt.add_edge(graph_components[component],art)
##            else:
##                ypages[node]=graph_components[component]
##    for art in articulations:
##        ypages[art] = art
def block_tree_find(bt,s,t):
    """
    Given a block tree bt and s and t articulation points sharing a
    biconnected components of the original graph, finds the
    given biconnected component in the block tree.
    """
    for G in bt.neighbors(s):
        if G.has_node(t):
            return G
    return None


def _other(edge,x):
    (u,v) = edge
    if u == x:
        return v
    elif v == x:
        return u
    else:
        raise Exception("Impossible to know which is the other edge since no edge coincides")

def bicomponents(G):
    """
    Finds the biconnected components and articulation points of G, being G a undirected graph.
    Runs in time O(|V|+|E|) and does not use recursion.
    Returns
      biconnected: a set containing the biconnected components, each expressed as a frozenset of nodes
      articulationPoints: a set containing the nodes that are articulation points.
    """
    no_of_nodes = G.number_of_nodes()
    block_tree = nx.Graph()
    low = {}
    num = {}
    path = []
    nextptr = {}
    edgestack = []
    articulation_points = set()
    found = {}
    components = set()
    for v in G:
        #print "A"
        low[v] = 0
        num[v] = 0
        found[v] = False
        nextptr[v] = 0

    for i in G:
        if low[i] != 0:
            continue
        path.append(i)
        counter = 1
        rootdfs = 0
        low[i] = counter
        num[i] = counter
        counter = counter + 1
        while len(path) > 0:
            act = path[-1]
        # print "act:",act,
            actnext = nextptr[act]
        # print "actnext",actnext,
            adjedges = G.edges(act)
        # print "adjedges",adjedges
            n = len(adjedges)
            if actnext < n:
        # Step down (maybe)
            # print "Stepping down"
                edge = adjedges[actnext]
            # print "Edge =",edge,
                nei = _other(edge,act)
            # print "=",(act,nei)
                if low[nei] == 0:
                    if act == i:
                        rootdfs = rootdfs + 1
                    # print "Increasing rootdfs to",rootdfs
                # print "We append",edge,"to",edgestack,
                    edgestack.append(edge)
                # print "and we get",edgestack
                    path.append(nei)
                    low[nei] = counter
                    num[nei] = counter
                    counter = counter + 1
                else:
                    # Update low value
                    low[act] = min(num[nei],low[act])
                nextptr[act] += 1
            else:
                # Step up
                path.pop()
                if not (len(path) == 0):
                    prev = path[-1]
                    # Update low value if needed
                    low[prev] = min(low[prev],low[act])
                    # Check for articulation point
                    if low[act] >= num[prev]:
                    # print "Should I add",prev," act is ",act, "i is",i
                        if not found[prev] and prev != i :
                        # print "YES!"
                            articulation_points.add(prev)
                            found[prev] = True
                        else:
                            pass
                        # print "NO"
                        #comp = nx.Graph()
                    # print edgestack
                    # print prev
                    # print "New component",
                        comp = set()
                        foundUV = False
                        while len(edgestack)>0:
                            (u,v) = edgestack.pop()
                            comp.add(u)
                            comp.add(v)
                            #comp.add_edge(u,v)
                        # print (u,v),
                            if (u == prev) or (v == prev):
                                break
                        #while len(edgestack) > 0:
                        #  (u,v) = edgestack[-1]
                        #  if u in comp.nodes() and v in comp.nodes():
                        #    edgestack.pop()
                        #    comp.add_edge(u,v)
                        #  else:
                        #    break
                    # print comp
                        components.add(frozenset(comp))
    # print "Should I add ",i ,"when rootdfs is",rootdfs
        if rootdfs >= 2:
        # print "YES!"
            articulation_points.add(i)
        else:
            pass
        # print "NO!"
    return (components,articulation_points)


def _show_block_list(bl,all_edges):
    s = 0
    for g in bl:
        edges = g
        s = s + len(edges)
        #print edges
    for (u,v) in all_edges:
        found = False
        for g in bl:
            if (u,v) in g or (v,u) in g:
                found = True
                break
        if not found:
            pass
            #print (u,v),"is not included"

    return s

def is_biconnected(g):
    nodes = g.nodes()
    for node in g:
        nodes.remove(node)
        h = g.subgraph(nodes)
        if not nx.algorithms.components.connected.is_connected(h):
            return False
        nodes.append(node)
    return True

def _check_all_biconnected(bt):
    for node in bt:
        if type(node)==nx.Graph:
            if not is_biconnected(node):
                return False
    return True

def _manually_find_articulation_points(G):
    nodes = G.nodes()
    articulation_points = set()
    num = nx.algorithms.components.connected.number_connected_components(G)
    for node in G:
        nodes.remove(node)
        H = G.subgraph(nodes)
        new_num = nx.algorithms.components.connected.number_connected_components(H)
        #print "for ",node,"the number is",
        if new_num > num:
            articulation_points.add(node)
        nodes.append(node)
    return articulation_points
