import networkx as nx
from itertools import combinations, permutations


class PhyNetwork:
    
    def __init__(self, G:nx.DiGraph): 
        self.digraph = nx.DiGraph()
        if (G is not None): 
            self.digraph.add_edges_from(G.edges)
            self.digraph.add_nodes_from(G.nodes)

        self._leaves = None
        self._reachables_by_node = None
        self._clusters_by_node = None
        self._hybrid_nodes = None
        self._root = None
        self._clustering_system = None
        self._overlap_graph = None
        self._hasse_diagram = None


    # returns list of reachable nodes of n
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
        
    # creates dict of reachable nodes for every node            
    def reachables_by_node(self):
        if self._reachables_by_node is not None: 
            return self._reachables_by_node

        reachables_by_node = dict()
        for node in self.digraph.nodes:
            reachables_by_node[node] = self.reachable_nodes(node)

        self._reachables_by_node = reachables_by_node
        return reachables_by_node

    # creates dict for clusters by node
    def clusters_by_node(self):
        if self._clusters_by_node is not None: 
            return self._clusters_by_node

        clusters_by_node = dict() 
        for node, reachables in self.reachables_by_node().items():
            reachable_leaves = set()
            for reachable in reachables: 
                if reachable in self.leaves():
                    reachable_leaves.add(reachable)
            clusters_by_node[node] = reachable_leaves

        self._clusters_by_node = clusters_by_node
        return clusters_by_node

    # returns the leaves of the digraph
    def leaves(self): 
        if self._leaves is not None: 
            return self._leaves

        leaves = set([leave for leave in self.digraph.nodes if self.digraph.out_degree(leave) == 0])

        self._leaves = leaves
        return leaves
    
    #returns the root of the digraph
    def root(self):
        if self._root is not None: 
            return self._root

        for node in self.digraph.nodes:
            if self.digraph.in_degree(node) == 0:
                self._root = node
                return node

    #returns all hybrid nodes of the digraph
    def hybrid_nodes(self): 
        if self._hybrid_nodes is not None:
            return self._hybrid_nodes

        hybrid_nodes = list()
        for node in self.digraph.nodes:
            if self.digraph.in_degree(node) > 1:
                hybrid_nodes.append(node)

        self._hybrid_nodes = hybrid_nodes
        return hybrid_nodes
        
    
    #creates and returns the overlap graph based on the clustering system of the digraph
    def overlap_graph(self):
        if self._overlap_graph is not None:
            return self._overlap_graph

        G = nx.Graph()
        for cluster_pair in combinations(self.clustering_system(),2): 
            intersection = cluster_pair[0].intersection(cluster_pair[1])
            #overlap condition 
            if intersection not in [set(), cluster_pair[0], cluster_pair[1]]: 
                # add edge between two clusters if they overlap
                G.add_edge(cluster_pair[0], cluster_pair[1], overlap = set(intersection))

        self._overlap_graph = G
        return G
    
    # creates the clustering system of the digraph    
    def clustering_system(self): 
        if self._clustering_system is not None:
            return self._clustering_system

        clustering_system = set()
        for cluster in self.clusters_by_node().values(): 
            clustering_system.add(frozenset(cluster))
        
        self._clustering_system = clustering_system
        return clustering_system
        

    # creates and returns the hasse diagram of the digraph based on the clusters
    def hasse_diagram(self):
        if self._hasse_diagram is not None:
            return self._hasse_diagram

        G = nx.DiGraph()
        
        for pair in permutations(self.clusters_by_node().values(),2):
            if pair[0].issubset(pair[1]) and pair[0] != pair [1]:
                skip = [pair[0], pair[1]]
                found_cluster_between = False
                for cluster in self.clusters_by_node().values():
                    if cluster in skip:
                        continue
                    elif (cluster.issubset(pair[1]) and pair[0].issubset(cluster)):
                        found_cluster_between = True
                if not found_cluster_between:
                    G.add_edge(frozenset(pair[1]),frozenset(pair[0]))
       
        self._hasse_diagram = G
        return G



