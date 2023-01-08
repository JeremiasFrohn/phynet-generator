import networkx as nx
from itertools import  combinations
import PhyNetwork as pn

# This Script uses results from Marc Hellmuth, David Schaller, and Peter F. Stadler. Clustering systems
# of phylogenetic networks, 2022


class PhyNetAnalyser: 
    def __init__(self, G:nx.DiGraph):
        self.tree_child = None
        self.shortcut_free = None
        self.phyNet = pn.PhyNetwork(G)
        self._is_pcc = None
        self._is_shortcut_free = None
        self._is_tree_child = None 
        self._level_k = None
        self._is_semi_regular = None
        self._is_prebinary = None
        
    #returns wether the digraph is a tree
    def is_tree(self):    
        if self.phyNet.hybrid_nodes():
                return False
        return True

    #returns wether the digraph is phylogenetic 
    def is_phylogenetic(self):
        if self.phyNet.root() is None: 
            return False
        for node in self.phyNet.digraph.nodes:
            if self.phyNet.digraph.out_degree(node) == 1 and self.phyNet.digraph.in_degree(node) <= 1:
                return False
        return True

    #returns the level of the digraph
    def level_k(self):
        if self._level_k is not None:
            return self._level_k

        max_hybrids = 0
        for biconnect_component in nx.biconnected_components(nx.Graph(self.phyNet.digraph.to_undirected())):
            current_hybrid_counter = 0
            for node in biconnect_component:
                if self.phyNet.digraph.in_degree(node) >= 2 and not set(self.phyNet.digraph.predecessors(node)).isdisjoint(biconnect_component):
                    #ignore highest hybrid 
                    #if any(n in set(G.predecessors(node)) for n in biconnect_component)
                    current_hybrid_counter = current_hybrid_counter + 1
            if current_hybrid_counter > max_hybrids:
                max_hybrids = current_hybrid_counter
        self._level_k = max_hybrids
        return max_hybrids


    #returns wether the digraph is tree child
    def is_tree_child(self):
        if self._is_tree_child is not None: 
            return self._is_tree_child

        tree_child = True
        
        for node in self.phyNet.digraph.nodes:
            # iterate over every node that is not a leave 
            if self.phyNet.digraph.out_degree(node) > 0:
                # set tree child false until at least one tree child could be found
                tree_child = False
                for child in self.phyNet.digraph.successors(node):
                    if self.phyNet.digraph.in_degree(child) == 1:
                        tree_child = True
                        break
                if not tree_child:
                    break
        self._is_tree_child = tree_child
        return tree_child

    #returns wether the digraph is pcc
    def is_pcc(self):
        # check if pcc property has already been calculated 
        if self._is_pcc is not None: 
            return self._is_pcc

        self._is_pcc = False
        reachables_by_node = self.phyNet.reachables_by_node()
        clusters_by_node = self.phyNet.clusters_by_node()
        for node1, node2 in combinations(self.phyNet.digraph.nodes,2):
            # check if nodes are comparable
            if (
                node1 in reachables_by_node[node2]
                or node2 in reachables_by_node[node1]
            ):
                # if the nodes are comparable then check for subset condition
                continue
            else: 
                # if nodes are not comparable check for subset
                if (
                    clusters_by_node[node1] <= clusters_by_node[node2] 
                    or clusters_by_node[node2] <= clusters_by_node[node1]
                ):
                    # if subset but not comparable graph is no pcc else continue 
                    self._is_pcc = False
                    return False
        self._is_pcc = True
        return True 
    
    #returns wether the digraph is shortcut free
    def is_shortcut_free(self):
        if self._is_shortcut_free is not None: 
            return self._is_shortcut_free

        self._is_shortcut_free = False

        reachables_by_node = self.phyNet.reachables_by_node()
        for node in self.phyNet.hybrid_nodes():
            #check parent nodes
            parent_nodes = self.phyNet.digraph.predecessors(node)
            for parent1, parent2 in combinations(parent_nodes,2): 
                if parent1 in reachables_by_node[parent2] or parent2 in reachables_by_node[parent1]:
                    return False
        self._is_shortcut_free = True                   
        return True 
        
    #returns wether the digraph is normal
    def is_normal(self):
        if self.is_shortcut_free() and self.is_tree_child():
            return True
        else:
            return False

    #returns wether the digraph is semi regular 
    def is_semi_regular(self):
        if self._is_semi_regular is not None: 
            return self._is_semi_regular

        if self.is_shortcut_free() and self.is_pcc():
            self._is_semi_regular = True
            return True
        else:
            self._is_semi_regular = False
            return False

    #returns wether the the digraph is regular
    def is_regular(self):
        #semiregular kein knoten mit outdegree 1 
        if self.is_semi_regular():
            for node in self.phyNet.digraph.nodes:
                if self.phyNet.digraph.out_degree(node) == 1:
                    return False
            return True
        return False

    #returns wether the digraph is separated
    def is_separated(self):
        for node in self.phyNet.hybrid_nodes():
            if not self.phyNet.digraph.out_degree(node) == 1:
                return False
        return True
    
    #returns wether the digraph is binary
    def is_binary(self):
        for node in self.phyNet.digraph.nodes:
            node_in = self.phyNet.digraph.in_degree(node)
            node_out = self.phyNet.digraph.out_degree(node)
            if (node_in,node_out) not in {(0,2),(1,0),(2,1),(1,2)}:
                return False
        return True

    #returns wether the digraph is cluster network
    def is_cluster_network(self):
        if (
            self.is_semi_regular()
            and self.is_separated()
            and self.is_phylogenetic()
        ):
            return True
        else:
            return False

    #returns wether the digraph is
    def is_galled_tree(self):

        if self.level_k() <= 1:
            for node in self.phyNet.hybrid_nodes():
                if self.phyNet.digraph.in_degree(node) != 2:
                    return False
            return True
        return False

   #returns wether the digraph is conventional
    def is_conventional(self):
        for leave in self.phyNet.leaves():
            if self.phyNet.digraph.in_degree(leave) > 1:
                return False

        # nx.biconnected_components yields sets of nodes, so they must contain more than two nodes to be non trivial 
        non_trivial_blocks = [bcc for bcc in nx.biconnected_components(nx.to_undirected(self.phyNet.digraph)) if len(bcc) > 2]
        for hybrid_node in self.phyNet.hybrid_nodes():
            # check if hybrid node is in unique  non_trivial_block
            unique = True
            # muss unique sein ! 
            for non_trivial_block in non_trivial_blocks: 
                if hybrid_node in non_trivial_block:
                    #continue with next hybrid node
                    if unique:
                        unique = False
                    else:
                        return False
                    
        return True

    #returns wether the digraph is quasi binary
    def is_quasi_binary(self):
        # in = 2 and out = 1 for every hybrid 
        for hybrid_node in self.phyNet.hybrid_nodes():
            if self.phyNet.digraph.out_degree(hybrid_node) != 1 or self.phyNet.digraph.in_degree(hybrid_node) != 2:
                return False

        # out (max B) = 2 for every non trivial block B 
        non_trivial_blocks = [bcc for bcc in nx.biconnected_components(nx.to_undirected(self.phyNet.digraph)) if len(bcc) > 2]
        for non_trivial_block in non_trivial_blocks:
            for node in non_trivial_block:
                if set(self.phyNet.digraph.predecessors(node)).isdisjoint(non_trivial_block):
                     if self.phyNet.digraph.out_degree(node) != 2:
                        return False
        return True

    #returns wether the clustering system is hierarchy
    def is_hierarchy(self):
        # if overlap graph is empty there are no overlaps, therefor, the clustering system is a hierarchy
        return nx.is_empty(self.phyNet.overlap_graph())

    # lemma 3.30
    # clustering system cl is closed if and only if A,B in CL AND A.intersection(B) not empty implies A.intersection(B) CL
    #returns wether the clustering system is closed
    def is_closed(self):
        # A and B in clustering system
        for cluster_A, cluster_B in combinations(self.phyNet.clustering_system(), 2):
            # if clusters A and B are not disjoint
            if not cluster_A.isdisjoint(cluster_B):
                # intersection of A and B is not in the clustering system
                if cluster_A.intersection(cluster_B) not in self.phyNet.clustering_system():
                    return False
        return True


    #Hasse Diagramm, for every leaf check if x,y is subset
    #for every x,y there is an inclusion minimal c with x,y subset of c
    #returns wether the clustering system is prebinary    
    def is_prebinary(self):
        if self._is_prebinary is not None:
            return self._is_prebinary

        G = self.phyNet.hasse_diagram()
        for leaf_pair in combinations(self.phyNet.leaves(),2):
            found = False
            for node in G.nodes:
                # length
                leaf_pair_set = set(leaf_pair)
                if leaf_pair_set.issubset(node):
                    if not [child for child in G.successors(node) if leaf_pair_set.issubset(child)]:
                        if found:
                            self._is_prebinary = False
                            return False
                        found = True
            if not found:
                self._is_prebinary = False
                return False   
        self._is_prebinary = True  
        return True
                        

    #returns wether the clustering system binary
    def is_binary_cl(self):
        if not self.is_prebinary():
            return False

        G = self.phyNet.hasse_diagram()

        for node in G:
            #skip singletons
            if len(node) == 1:
                continue
            #clusters with two leaves x,y are inclusion minimal for x,y
            if len(node) == 2:
                continue
            found = False
            for leaf_pair in combinations(node,2):
                # if children do not contain leaf_pair then 
                if not [child for child in G.successors(node) if set(leaf_pair).issubset(child)]:
                    found = True
                    break
            if not found:
                return False
        return True

    #returns wether the clustering system is weak hierarchy
    def is_weak_hierarchy(self):
        for triple in combinations(self.phyNet.clustering_system(), 3): 
            intersectionC1_C2 = frozenset(triple[0].intersection(triple[1]))
            intersectionC1_C3 = frozenset(triple[0].intersection(triple[2]))
            intersectionC2_C3 = frozenset(triple[1].intersection(triple[2]))
            if (intersectionC1_C2.intersection(triple[2])
                in set([intersectionC1_C2, intersectionC2_C3, intersectionC1_C3,frozenset()])):
                continue
            else:
                return False
        return True

    #returns wether the network is L
    def is_L(self):
        overlap_graph = self.phyNet.overlap_graph()
        for node in overlap_graph.nodes:
            neighbors = [neighbor for neighbor in overlap_graph.neighbors(node)]
            if len(neighbors) >=2:
                # edge data cluster of node and first neighbor
                overlap = overlap_graph[node][neighbors[0]]['overlap']
                for i in range(1,len(neighbors)):
                    if overlap == overlap_graph[node][neighbors[i]]['overlap']:
                        continue
                    else:
                        return False          
        return True

    #returns wether the network is N30
    def is_N30(self):
        # nx.triangles returns dictionary with count of all triangles for every node in G, a triangle is counted three times
        for value in nx.triangles(self.phyNet.overlap_graph()).values():
            if value != 0:
                return False
        return True

    # every cluster overlaps with at most one other cluster 
    #returns wether the clustering system is paired hierarchy   
    def is_paired_hierarchy(self):
        overlap_graph = self.phyNet.overlap_graph()
        for node in overlap_graph:
            if overlap_graph.degree(node) > 1:
                return False
        return True

    #returns wether the network is 2inc
    def is_2_inc(self): 
        G = self.phyNet.hasse_diagram()
        # in and out at max 2 
        for node in G.nodes:
            if G.out_degree(node) <= 2 and G.in_degree(node) <= 2:
                continue
            else:
                return False 
        return True 
        



