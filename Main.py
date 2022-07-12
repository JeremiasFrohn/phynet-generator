from itertools import pairwise
import networkx as nx
import matplotlib.pyplot as plt
import asymmetree.treeevolve as te
from asymmetree.tools.PhyloTreeTools import to_newick
from tralda.datastructures.Tree import Tree, TreeNode
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
import random
import numpy
import vistools as vt
import Simulator as sim
import PhyNetAnalyser as pna
import timeit

s1 = sim.Simulator()
s2 = sim.Simulator()
G = s1.simulate(5, 1, 2, 1, 1, 0)
G2 = s2.simulate_by_tree(3)
s1.create_cluster_dict()
print("cluster dict")
print(s1.cluster_dict)
print("successors of G")
print(s1.node_dict)
print(s2.node_dict)

pos = vt.topo_pos(G)
pos2 = vt.topo_pos(G2)

print("reachable nodes from 1")
node_list = s1.reachable_nodes(1)
print(node_list)

G3 = nx.Graph(s1.digraph.to_undirected())

print("biconnected components")
for biconnected_component in nx.biconnected_components(G3):
    print("bc:")
    print(biconnected_component)

analyser = pna.PhyNetAnalyser(G)
print("is shortcut free:")
print(str(analyser.is_shortcut_free(G)))

print("is pcc:")
print(str(analyser.is_pcc(G)))

print("return leave function ")
print(s1.return_leaves())

print("predecessors of 2 ")
print(set(s1.digraph.predecessors(2)))


liste = [1,2,3,4,5]
print("pairwise list iteration")
for e1, e2 in pairwise(liste):
    print(str([e1,e2]))


plt.subplot(221)
nx.draw_networkx(
    G,
    pos=pos,
    arrows=True,
    with_labels=True,
)
plt.subplot(222)
nx.draw_networkx(
    G2,
    pos=pos2,
    arrows=True,
    with_labels=True,
)
plt.show()
