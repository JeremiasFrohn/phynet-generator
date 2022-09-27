
from ast import List
from asyncio.windows_events import NULL
import asymmetree.treeevolve as te
import networkx as nx
import PhyNetwork as pn
import random
import itertools

class Simulator:
    def __init__(self):
        self.phyNet = pn.PhyNetwork(nx.DiGraph())


    # weight, oder h raus und m type und n type sind die summen 
    def simulate(self, n, s, h, mtype_weight, ntype_weight):
        # initalize Graph with root and one specialization
        self.phyNet.digraph.add_edges_from([(1, 2), (1, 3)])
        # leaves that will be updated whenever the leaves change
        leaves = list()
        leaves.extend([2,3])
        # node and leave counter
        node_counter = 4
        leave_counter = 2
        # sum of all hybridization weights for
        weight = ntype_weight + mtype_weight 
        # specialization and hybridization process
        while len(leaves) < n:
            # specialization
            rdm = random.random()
            rdm2 = random.random()
            if rdm <= s / (h + s):
                x = random.choice(leaves)
                self.phyNet.digraph.add_edges_from([(x, node_counter), (x, node_counter + 1)])
                leaves.remove(x)
                leaves.extend([node_counter, node_counter + 1])
                leave_counter += 2
                node_counter = node_counter + 2
            # n_type hybridization
            elif rdm2 <= ntype_weight / weight:
                xy = random.sample(leaves, 2)
                leaves.remove(xy[0])
                leaves.remove(xy[1])
                leaves.extend([node_counter, node_counter + 1])
                self.phyNet.digraph.add_edges_from(
                    [(xy[0], node_counter), (xy[1], node_counter + 1), (xy[0], xy[1])]
                )
                node_counter = node_counter + 2
            # m_type hybridization
            else:
                xy = random.sample(leaves, 2)
                leaves.remove(xy[0])
                leaves.remove(xy[1])
                leaves.extend([node_counter, node_counter + 1, node_counter + 2])

                self.phyNet.digraph.add_edges_from(
                    [
                        (xy[0], node_counter),
                        (xy[1], node_counter + 1),
                        (xy[0], node_counter + 2),
                        (xy[1], node_counter + 2),
                    ]
                )
                node_counter = node_counter + 3
                leave_counter += 1
        return self.phyNet.digraph

    def simulate_by_tree(self, n, k):
        tree = te.species_tree_n(n, planted=False)
        #node dict erstellen udn updaten zum hinzufügen von kanten
        for edge in tree.edges():
            self.phyNet.digraph.add_edge(edge[0].label, edge[1].label)
        
        for i in range(k):
            self.add_edge_random()
        
        return self.phyNet.digraph

    def add_edge_random(self, subgraph_edges = None):
        # take 2 random edges e1, e2 and split them via two nodes v1, v2
        edges = list()
        if subgraph_edges is not None:
            edges = list(subgraph_edges)
        else:
            edges = [e for e in self.phyNet.digraph.edges]
        # sample
        e1,e2 = random.sample(edges,2)
        max_node= max([node for node in self.phyNet.digraph.nodes()])
        v1 = max_node + 1
        v2 = v1 + 1

        self.phyNet.digraph.remove_edges_from([e1, e2]) 
        self.phyNet.digraph.add_edges_from([(e1[0], v1), (v1, e1[1])])
        self.phyNet.digraph.add_edges_from([(e2[0], v2), (v2, e2[1])])

        # if v1 is above v2 add edge from v1 to v2 else add edge from v2 to v1
        if nx.has_path(self.phyNet.digraph, v2, v1):
            self.phyNet.digraph.add_edge(v2, v1)
        else:
            self.phyNet.digraph.add_edge(v1, v2)

    def simulate_bcc(self, n, contraction_prob, binary, level_k = 1):
        #catch error
        if level_k <= 0:
            print("cant simulate level k smaller than 1")
            return
        
        tree = te.species_tree_n(n, planted = False, contraction_probability = contraction_prob)
        for edge in tree.edges():
            self.phyNet.digraph.add_edge(edge[0].label, edge[1].label)

        #list all non binary nodes
        nonbinary_nodes = [node for node in self.phyNet.digraph.nodes if self.phyNet.digraph.out_degree(node) > 2]
       

        for non_binary_node in nonbinary_nodes:
            out_degree = self.phyNet.digraph.out_degree(non_binary_node)
            node_counter = max([node for node in self.phyNet.digraph.nodes()]) + 1
            #remember children
            children = [child for child in self.phyNet.digraph.successors(non_binary_node)]
            #remove edges from non binary node to children
            for child in children:
                self.phyNet.digraph.remove_edge(non_binary_node, child)
            
            hybrid_node = node_counter + out_degree - 1
            nodes_to_add = [node for node in range(node_counter, node_counter+out_degree-1)]
            #slice
            nodes_for_level_k = [node for node in nodes_to_add]
            nodes_for_children = [node for node in nodes_to_add]
            nodes_for_children.append(hybrid_node)
            lowest_nodes = list()
            lowest_nodes.append(non_binary_node)
            
            #set outdegree for root of the biconnected component
            if binary:
                out = 2
            else: 
                out = out_degree 

            while nodes_to_add:
                # chose nodes for new edge
                random_lowest_node = random.choice(lowest_nodes)
                random_node_to_add = random.choice(nodes_to_add)
                # add new edge
                self.phyNet.digraph.add_edge(random_lowest_node, random_node_to_add)
                #update lowest nodes
                if random_lowest_node is not non_binary_node or self.phyNet.digraph.out_degree(non_binary_node) == out:
                    lowest_nodes.remove(random_lowest_node)  
                lowest_nodes.append(random_node_to_add)
                nodes_to_add.remove(random_node_to_add)
            

            for lowest_node in lowest_nodes:
                self.phyNet.digraph.add_edge(lowest_node, hybrid_node)
            
            for node in nodes_for_children:
                self.phyNet.digraph.add_edge(node, children.pop())

            # level_k 
            # kanten auflösen
            for i in range(level_k - 1):
                self.add_edge_random(nx.edge_bfs(self.phyNet.digraph, non_binary_node))
                    


        return self.phyNet.digraph

    def simulate_non_binary(self, n, contraction_prob):
        #non binary tree to digraph
        print("create phyNet by Tree")
        tree = te.species_tree_n(n, planted = False, contraction_probability = contraction_prob)
        for edge in tree.edges():
            self.phyNet.digraph.add_edge(edge[0].label, edge[1].label)

        # vorher liste der nb knoten machen 
        nonbinary_nodes = [node for node in self.phyNet.digraph.nodes() if self.phyNet.digraph.out_degree(node) > 2]
        node_counter = max([node for node in self.phyNet.digraph.nodes()]) + 1
        for node in nonbinary_nodes:
            # for node with outdegree > 2 create new tree 
            out_degree = self.phyNet.digraph.out_degree(node)
            
            print("Node with outdegree > 2 found")
            print(out_degree)
            # remember the children of non binary node
            children = [child for child in self.phyNet.digraph.successors(node)]
            # remove edge to children
            for child in children:
                self.phyNet.digraph.remove_edge(node, child)               
            #simulate tree that gets connected to the binary node
            tree1 = te.species_tree_n( (out_degree - 1) , planted=False)
            leaves_tree1 = list(tree1.leaves())
            #count nodes for the labels
            
            # update labels
            for tree_node in tree1.preorder():
                if tree_node is tree1.root:
                    #set root lable to nb node
                    tree1.root.label = node
                else:
                    tree_node.label = node_counter
                    node_counter = node_counter + 1
            # add edges to digraph
            for edge in tree1.edges():
                self.phyNet.digraph.add_edge(edge[0].label, edge[1].label)
            
            for leave in leaves_tree1:
                self.phyNet.digraph.add_edge(leave.label, children.pop())

            #simulate new tree with offset 
            tree2 = te.species_tree_n((out_degree - 1), planted = False, contraction_probability = contraction_prob)
            
            for tree_node in tree2.postorder():
                if tree_node.is_leaf():
                    tree_node.label = leaves_tree1.pop().label
                else:
                    tree_node.label = node_counter
                    node_counter = node_counter + 1
                    
            for edge in tree2.edges():
                self.phyNet.digraph.add_edge(edge[1].label, edge[0].label)

            self.phyNet.digraph.add_edge(tree2.root.label, children.pop())
                

        return self.phyNet.digraph



        #return tuple consisting of edges of each split




    




    