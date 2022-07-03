import asymmetree.treeevolve as te
import networkx as nx
import random

class Simulator:
    def __init__(self):
        self.digraph = nx.DiGraph()
        self.leaves = []
        self.node_dict = {}
        self.cluster_dict = {}

    def simulate(self, n, s, h, mtype_rate, ntype_rate, ytype_rate):
        # initalize Graph with root and one specialization
        self.digraph.add_edges_from([(1, 2), (1, 3)])
        # keep track of reachable nodes
        self.node_dict[1] = set([2, 3])
        # leaves that will be updated whenever the leaves change
        self.leaves.extend([2, 3])
        # node and leave counter
        node_counter = 4
        leave_counter = 2
        # sum of all hybridization rates for
        rate_sum = ntype_rate + mtype_rate + ytype_rate
        # specialization and hybridization process
        while len(self.leaves) < n:
            # specialization
            rdm = random.random()
            rdm2 = random.random()
            if rdm <= s / (h + s):
                print("added specialisation")
                x = random.choice(self.leaves)
                self.digraph.add_edges_from([(x, node_counter), (x, node_counter + 1)])
                self.update_node_dict(node_counter)
                self.update_node_dict(node_counter + 1)
                self.leaves.remove(x)
                self.leaves.extend([node_counter, node_counter + 1])
                leave_counter += 2
                node_counter = node_counter + 2
            # n_type hybridization
            elif rdm2 <= ntype_rate / rate_sum:
                print("added ntype")
                xy = random.sample(self.leaves, 2)
                self.leaves.remove(xy[0])
                self.leaves.remove(xy[1])
                self.leaves.extend([node_counter, node_counter + 1])
                self.digraph.add_edges_from(
                    [(xy[0], node_counter), (xy[1], node_counter + 1), (xy[0], xy[1])]
                )
                self.update_node_dict(node_counter)
                self.update_node_dict(node_counter + 1)
                node_counter = node_counter + 2
            # m_type hybridization
            else:
                print("added mtype")
                xy = random.sample(self.leaves, 2)
                self.leaves.remove(xy[0])
                self.leaves.remove(xy[1])
                self.leaves.extend([node_counter, node_counter + 1, node_counter + 2])

                self.digraph.add_edges_from(
                    [
                        (xy[0], node_counter),
                        (xy[1], node_counter + 1),
                        (xy[0], node_counter + 2),
                        (xy[1], node_counter + 2),
                    ]
                )
                self.update_node_dict(node_counter)
                self.update_node_dict(node_counter + 1)
                self.update_node_dict(node_counter + 2)
                node_counter = node_counter + 3
                leave_counter += 1
        print("Edgelist sim Graph")
        print(self.digraph.edges)
        return self.digraph

    def simulate_by_tree(self, n):
        tree = te.simulate_species_tree(n, planted=False)
        for edge in tree.edges():
            self.digraph.add_edge(edge[0].label, edge[1].label)
        
        for i in range(n):
            self.add_edge_random()
        print("Edgelist sim by tree Graph")
        print(self.digraph.edges)
        self.create_node_dict()
        return self.digraph

    def add_edge_random(self):
        # take 2 random edges e1, e2 and split them via two nodes v1, v2
        edges = [e for e in self.digraph.edges]
        # sample
        e1 = random.choice(edges)
        e2 = random.choice(edges)
        v1 = len(self.digraph)
        v2 = v1 + 1

        self.digraph.remove_edges_from([e1, e2])
        self.digraph.add_edges_from([(e1[0], v1), (v1, e1[1])])
        self.digraph.add_edges_from([(e2[0], v2), (v2, e2[1])])

        # if v1 is above v2 add edge from v1 to v2 else add edge from v2 to v1
        if nx.has_path(self.digraph, v2, v1):
            self.digraph.add_edge(v2, v1)
        else:
            self.digraph.add_edge(v1, v2)

    #TO-DO benchmark vs reachable nodes with outdegree
    def reachable_nodes_leaves(self, n):
        leaves = set(self.leaves)
        # list of nodes that need to be checked for succcessors
        nodes_to_check = [successor for successor in self.digraph.successors(n)]
        # list to return
        reachable_nodes = set(nodes_to_check)
        while len(nodes_to_check) != 0:
            # if all nodes in node to check are leaves the while loop shall stop
            for node in nodes_to_check:
                # check node for successors
                if node not in leaves:  
                    # if node has successors then add the successors to the reachable nodes and the list nodes to check
                    for successor in self.digraph.successors(node):
                        #if successor is already in reachable nodes continue with next successor 
                        if successor in reachable_nodes:
                            continue
                        else:
                            #add the successor to the nodes to check and the reachable nodes
                            nodes_to_check.append(successor)
                            reachable_nodes.add(successor)       
                # remove node because it has been checked for successors
                nodes_to_check.remove(node)
        return reachable_nodes


    def reachable_nodes(self, n):
        # list of nodes that need to be checked for succcessors
        nodes_to_check = [successor for successor in self.digraph.successors(n)]
        # list to return
        reachable_nodes = set(nodes_to_check)
        while len(nodes_to_check) != 0:
            # if all nodes in node to check are leaves the while loop shall stop
            for node in nodes_to_check:
                # check node for successors
                if self.digraph.out_degree(node) != 0:  # eventuell ersetzen mit leaves set zur schnelleren abfrage
                    # if node has successors then add the successors to the reachable nodes and the list nodes to check
                    for successor in self.digraph.successors(node):
                        #if successor is already in reachable nodes continue with next successor 
                        if successor in reachable_nodes:
                            continue
                        else:
                            #add the successor to the nodes to check and the reachable nodes
                            nodes_to_check.append(successor)
                            reachable_nodes.add(successor)       
                # remove node because it has been checked for successors
                nodes_to_check.remove(node)
        return reachable_nodes


    def update_node_dict(self, n):
        """Updates the node dictory whenever there is a new specialisation/hybridization during the simulation"""
        self.node_dict[n] = set()
        nodes_to_check = [predecessor for predecessor in self.digraph.predecessors(n)]
        for node in nodes_to_check:
            if node in self.node_dict:
                self.node_dict[node].add(n)
            else:
                self.node_dict[node] = set([n])
        
        while len(nodes_to_check) != 0:
            for node in nodes_to_check:
                for predecessor in self.digraph.predecessors(node):
                    nodes_to_check.append(predecessor)
                    if predecessor in self.node_dict:
                        self.node_dict[predecessor].add(n)
                    elif self.digraph.out_degree(predecessor) == 0:
                            self.node_dict[predecessor] = set()
                    else:
                        self.node_dict[predecessor] = set([n])
                nodes_to_check.remove(node)
        
    def create_node_dict(self):
        for node in self.digraph.nodes:
            self.node_dict[node] = self.reachable_nodes(node)

    
    def create_cluster_dict(self):
        for node, reachables in self.node_dict.items():
            reachable_leaves = set()
            for reachable in reachables: 
                if self.digraph.out_degree(reachable) == 0:
                    reachable_leaves.add(reachable)
            self.cluster_dict[node] = reachable_leaves

    
    def return_leaves(self):
        return set([leave for leave in self.digraph.nodes if self.digraph.out_degree(leave) == 0])
        

    def level_k(self):
        max_hybrids = 0
        G = nx.Graph(self.digraph.to_undirected())
        for biconnect_component in nx.biconnected_components(G):
            current_hybrid_counter = 0
            for node in biconnect_component:
                if self.digraph.in_degree(node) >= 2:
                    current_hybrid_counter = current_hybrid_counter + 1
            if current_hybrid_counter > max_hybrids:
                max_hybrids = current_hybrid_counter
        return max_hybrids
        



    