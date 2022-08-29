from copy import copy
import networkx as nx
from itertools import combinations


class PhyNetwork:
    
    def __init__(self, G:nx.DiGraph): 
        self.digraph = nx.DiGraph()
        if (G is not None): 
            self.digraph.add_edges_from(G.edges)
            self.digraph.add_nodes_from(G.nodes)

        self.leaves = self.create_leave_set()
        self.reachables_by_node = self.create_node_dict()
        self.clusters_by_node = self.create_cluster_dict()
        self.root = self.find_root()
        self.clustering_system = self.create_clustering_system()
        #clustering system


    def reachable_nodes(self, n):
        # list of nodes that need to be checked for succcessors
        nodes_to_check = [successor for successor in self.digraph.successors(n)]
        nodes_to_check_copy = list(nodes_to_check)
        # list to return
        reachable_nodes = set(nodes_to_check)
        reachable_nodes.add(n)
        while nodes_to_check:
            # if all nodes in node to check are leaves the while loop shall stop
            for node in nodes_to_check:  
                    # if node has successors then add the successors to the reachable nodes and the list nodes to check
                for successor in self.digraph.successors(node):
                    #if successor is already in reachable nodes continue with next successor 
                    if successor in reachable_nodes:
                        continue
                    else:
                        #add the successor to the nodes to check and the reachable nodes
                        nodes_to_check_copy.append(successor)
                        reachable_nodes.add(successor)       
                # remove node because it has been checked for successors
                nodes_to_check_copy.remove(node)
            #overwrite nodes_to_check with the updated list
            nodes_to_check = list(nodes_to_check_copy)
        return reachable_nodes
        
                
    def create_node_dict(self):
        reachables_by_node = dict()
        for node in self.digraph.nodes:
            reachables_by_node[node] = self.reachable_nodes(node)
        return reachables_by_node

    
    def create_cluster_dict(self):
        clusters_by_node = dict() 
        for node, reachables in self.reachables_by_node.items():
            reachable_leaves = set()
            for reachable in reachables: 
                if self.digraph.out_degree(reachable) == 0:
                    reachable_leaves.add(reachable)
            clusters_by_node[node] = reachable_leaves
        return clusters_by_node

    def create_leave_set(self): 
        return set([leave for leave in self.digraph.nodes if self.digraph.out_degree(leave) == 0])


    def leaves(self):
        for node in self.digraph.nodes:
            if self.digraph.out_degree(node) == 0:
                yield node


    def cluster(self, n):
        #generate all leaves reachable from node n 
        for node in self.reachable_nodes(n):
            if self.digraph.out_degree(node) == 0:
                yield node
    

    def find_root(self):
        for node in self.digraph.nodes:
            if self.digraph.in_degree(node) == 0:
                return node
        
        
    def overlap_graph(self):
        if self.clusters_by_node is None: 
            self.create_cluster_dict()

        G = nx.Graph()
        for cluster_pair in combinations(self.clusters_by_node.values(),2 ) : 
            intersection = cluster_pair[0].intersection(cluster_pair[1])
            #overlap condition 
            if intersection not in [set(), cluster_pair[0], cluster_pair[1]]: 
                # add edge between two clusters if they overlap
                G.add_edge(frozenset(cluster_pair[0]), frozenset(cluster_pair[1]), overlap = intersection)

        print("overlap graph")
        
        print(G.edges)
        nx.get_edge_attributes(G, "overlap")
        return G
    
        
    def create_clustering_system(self): 
        if self.clusters_by_node is None: 
            self.clusters_by_node = self.create_cluster_dict
        
        clustering_system = set()
        for cluster in self.clusters_by_node.values(): 
            clustering_system.add(frozenset(cluster))
        
        return clustering_system
        
