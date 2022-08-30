from itertools import pairwise

import networkx as nx
import matplotlib.pyplot as plt
import asymmetree.treeevolve as te
import vistools as vt
import Simulator as sim
import PhyNetwork
import PhyNetAnalyser as pna
import asymmetree


s1 = sim.Simulator()
s2 = sim.Simulator()
s3 = sim.Simulator()
s4 = sim.Simulator()
G = s1.simulate(5, 1, 2, 1, 1)
G2 = s2.simulate_by_tree(3,2)
G3 = s3.simulate_non_binary(5, 0.8)
G4 = s4.simulate_bcc(5, 0.8, level_k=1, binary=False)
graph = nx.DiGraph()


# G4 = nx.DiGraph()
# tree = te.species_tree_n(8, planted = False, contraction_proportion = 0.9)
# for edge in tree.edges():
#     G4.add_edge(edge[0].label, edge[1].label)
# print("edges tree")
# print(G4.edges())
# print("cluster dict")
# print(s1.cluster_dict)
# print("successors of G")
# print(s1.node_dict)
# print(s2.node_dict)

#pos = vt.topo_pos(G)
#pos2 = vt.topo_pos(G2)
pos3 = vt.topo_pos(G3)
pos4 = vt.topo_pos(G)


print("G4")
print(G)
analyseG4 = pna.PhyNetAnalyser(G)
print(G)
print("Graph is closed:")
print(analyseG4.is_closed())
print("Graph is (L):")
print(analyseG4.is_L())
#overlap_G = analyseG4.phyNet.create_overlap_graph()
#print("hasse diagramm ")
#G5 = nx.DiGraph(analyseG4.hasse_diagramm())
#print(G5.edges)
#pos5 = vt.topo_pos(G5)


# analyser = pna.PhyNetAnalyser(G)
# print("is shortcut free:")
# print(str(analyser.is_shortcut_free(G)))

# print("is pcc:")
# print(str(analyser.is_pcc(G)))

# print("return leave function ")
# print(s1.return_leaves())

# print("predecessors of 2 ")
# print(set(s1.digraph.predecessors(2)))


# liste = [1,2,3,4,5]
# print("pairwise list iteration")
# for e1, e2 in pairwise(liste):
#     print(str([e1,e2]))

print("edges from g3 simulate non binary")
print(s4.phyNet.digraph.edges())
# plt.subplot(221)
# nx.draw_networkx(
#     G,
#     pos=pos,
#     arrows=True,
#     with_labels=True,
# )
# plt.subplot(222)
# nx.draw_networkx(
#     G3,
#     pos=pos3,
#     arrows=True,
#     with_labels=True,
# )


color_map = []

for node in G: 
    if G.in_degree(node) >= 2:
        color_map.append('red')
    else:
        color_map.append('blue')

#plt.figure(1)

spec_G = nx.DiGraph()
spec_G.add_edge(1,2)
spec_G.add_edge(1,3)
spec_pos = vt.topo_pos(spec_G)

nx.draw_networkx(
    G,
    pos = pos4,
    arrows=True,
    with_labels=True,
    node_color = color_map,
)

# nx.draw_networkx(
#     spec_G,
#     pos = spec_pos,
#     arrows=True,
#     with_labels=True,
# )

plt.savefig("Graph.png", format="png")
# plt.figure(2)
# pos = nx.spring_layout(overlap_G)

# edge_labels = nx.get_edge_attributes(overlap_G,'overlap')
# nx.draw_networkx_edge_labels(overlap_G, pos, edge_labels)

# nx.draw_networkx(
#     overlap_G,
#     pos
# )

plt.axis('off')
plt.show()

def simulate_multiple(n, nodes, s, h, mtype, ntype): 
    s = sim.Simulator()
    G = nx.DiGraph()
    G = s.simulate()
    