import Simulator
import PhyNetAnalyser
import matplotlib.pyplot as plt
import vistools as vt
import networkx as nx
import csv

def not_binary():
    no_binary_not_found = True
    while no_binary_not_found:
        simulator = Simulator.Simulator()
        G = simulator.simulate_by_tree(20,5)
        ana_G = PhyNetAnalyser.PhyNetAnalyser(G)

        if not ana_G.is_binary() :
            no_binary_not_found = False

            color_map = []

            for node in G: 
                if node in ana_G.phyNet.hybrid_nodes():
                    color_map.append('red')
                else:
                    color_map.append('blue')
            print(G.edges)
            nx.draw_networkx(
                G,
                pos = vt.topo_pos(G),
                arrows=True,
                with_labels=True,
                node_color = color_map
            )
            plt.show()

not_binary()