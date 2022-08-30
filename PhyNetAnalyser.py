from xml.dom.minidom import Element
import networkx as nx
from itertools import chain, combinations, permutations
import PhyNetwork as pn


class PhyNetAnalyser: 
    def __init__(self, G:nx.DiGraph):
        #self.tree_child = None
        #self.shortcut_free = None
        self.phyNet = pn.PhyNetwork(G)
        self.is_pcc = None
        self.is_shortcut_free = None
        self.is_tree_child = None 
        self.level_k = None
        self.is_semi_regular = None
        
    def is_tree(self):    
        for node in self.phyNet.digraph.nodes:
            if self.phyNet.digraph.in_degree(node) > 1:
                return True
        return False

    def is_phylogenetic(self):
        if self.phyNet.root is None: 
            return False
        for node in self.phyNet.digraph.nodes:
            if self.phyNet.digraph.out_degree(node) == 1 and self.phyNet.digraph.in_degree(node) <= 1:
                return False
        return True


    def level_k(self):
        if self.level_k is not None:
            return self.level_k

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
        self.level_k = max_hybrids
        return max_hybrids



    def is_tree_child(self):
        if self.is_tree_child is not None: 
            return self.is_tree_child

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
        self.is_tree_child = tree_child
        return tree_child


    def is_pcc(self):
        # check if pcc property has already been calculated 
        if self.is_pcc is not None: 
            return self.is_pcc

        self.is_pcc = False
        for pair in permutations(self.phyNet.digraph.nodes):
            # check if nodes are comparable
            if (
                pair[0] in self.phyNet.reachables_by_node[pair[1]]
                or pair[1] in self.phyNet.reachables_by_node[pair[0]]
            ):
                # if the nodes are comparable then check for subset condition
                if (
                    self.phyNet.clusters_by_node[pair[0]] <= self.phyNet.clusters_by_node[pair[1]]  
                    or self.phyNet.clusters_by_node[pair[1]] <= self.phyNet.clusters_by_node[pair[0]]
                ):
                    # if comparable and subset continue
                    continue
                else:
                    # if comparable but not subset graph is no pcc
                    print(self.phyNet.reachables_by_node)
                    print("comparable nodes where subset condition fails")
                    print([pair[0],pair[1]])
                    self.is_pcc = False
                    return False
            else: 
                # if nodes are not comparable check for subset
                if (
                    self.phyNet.clusters_by_node[pair[0]] <= self.phyNet.clusters_by_node[pair[1]] 
                    or self.phyNet.clusters_by_node[pair[1]] <= self.phyNet.clusters_by_node[pair[0]]
                ):
                    # if subset but not comparable graph is no pcc else continue 
                    print(self.phyNet.reachables_by_node)
                    print("non comparable nodes where subset condition is true")
                    print([pair[0],pair[1]])
                    self.is_pcc = False
                    return False
        self.is_pcc = True
        return True 
    
    def is_shortcut_free(self):
        if self.is_shortcut_free is not None: 
            return self.is_shortcut_free

        self.is_shortcut_free = False
        for node in self.phyNet.digraph.nodes:
            # iterate over hybrid nodes only
            if self.phyNet.digraph.in_degree(node) >= 2:
                predecessors = self.phyNet.digraph.predecessors(node)
                for predecessor in predecessors:
                    for child in self.phyNet.digraph.successors(predecessor):
                    # check for every predecessor of the hybrid node if there is a child that is not the hybrid node itself
                    # but there is also a path from the child of the predecessor to the hybrid node
                        if child is not node and nx.has_path(self.phyNet.digraph, child, node) :
                            return False     
        self.is_shortcut_free = True                   
        return True 
        

    def is_normal(self,G):
        if self.is_shortcut_free() and self.is_tree_child():
            return True
        else:
            return False


    def is_semi_regular(self):
        if self.is_semi_regular is not None: 
            return self.is_semi_regular

        if self.is_shortcut_free() and self.is_pcc():
            self.is_semi_regular = True
            return True
        else:
            self.is_semi_regular = False
            return False
    
    def is_regular(self):
        #semiregular kein knoten mit outdegree 1 
        if self.is_semi_regular():
            for node in self.phyNet.digraph.nodes:
                if self.phyNet.digraph.out_degree(node) == 1:
                    return False
        return True


    def is_separated(self):
        # if rausziehen
        for node in [node for node in self.phyNet.digraph.nodes if self.phyNet.digraph.in_degree(node) >= 2] :
            if not self.phyNet.digraph.out_degree(node) == 1:
                return False
        return True
    
    def is_binary(self):
        # hybrid mit in 3 
        for node in [node for node in self.phyNet.digraph.nodes if self.phyNet.digraph.in_degree(node) <= 1]:
            if not self.phyNet.digraph.out_degree(node) in [0,2]:
                return False
        if self.is_separated():
            return True
        else:
            return False

    # Proposition 4.23 
    def is_cluster_network(self):
        if (
            self.is_semi_regular()
            and self.is_separated()
            and self.is_phylogenetic()
        ):
            return True
        else:
            return False

    #Corollary 8.5
    def is_galled_tree(self):
        if self.level_k() == 1:
            for node in self.phyNet.digraph.nodes:
                if self.phyNet.digraph.in_degree(node) != 2:
                    return False
        return True


    def is_conventional(self):
        for leave in self.phyNet.leaves():
            if self.phyNet.digraph.in_degree(leave) != 1:
                return False
        is_conventional = False
        # nx.biconnected_components yields sets of nodes, so they must contain more than two nodes to be non trivial 
        non_trivial_blocks = [bcc for bcc in nx.biconnected_components(nx.to_undirected(self.phyNet.digraph)) if len(bcc) >= 2]
        for hybrid_node in self.phyNet.hybrid_nodes:
            # check if hybrid node is in non_trivial_block
            is_conventional = False
            for non_trivial_block in non_trivial_blocks: 
                if hybrid_node in non_trivial_block:
                    #continue with next hybrid node
                    is_conventional = True
                    break
        return is_conventional


    def is_quasi_binary(self):
        # in = 2 and out = 1 for every hybrid 
        for hybrid_node in self.phyNet.hybrid_nodes:
            if self.phyNet.digraph.out_degree(hybrid_node) != 1:
                return False

        # out (max B) = 2 for every non trivial block B 
        non_trivial_blocks = [bcc for bcc in nx.biconnected_components(nx.to_undirected(self.phyNet.digraph)) if len(bcc) >= 2]
        for non_trivial_block in non_trivial_blocks:
            for node in non_trivial_blocks:
                if set(self.phyNet.digraph.predecessors(node)).isdisjoint(non_trivial_block):
                     if self.phyNet.digraph.out_degree(node) != 2:
                        return False
        return True





    # lemma 3.30
    # clustering system cl is closed if and only if A,B in cl and a intersection b not empty implies a intersection b in cl
    def is_closed(self):
        clustering_system = self.phyNet.clustering_system
        clustering_system_with_empty_set = clustering_system.union(set(set()))
        subsets_cs = self._powerset(clustering_system, allow_empty_set=False)
        for subset in subsets_cs: 
            if set.intersection(*[set(element) for element in subset]) not in clustering_system_with_empty_set:
                return False
        return True

    #overlap graph
    def is_L(self):
        clustering_system = set([frozenset(cluster) for cluster in self.phyNet.clusters_by_node.values()])
        for triple in permutations(clustering_system, 3):
            if (
                triple[0].intersection(triple[1]) not in {frozenset(), triple[0], triple[1]} 
                and triple[0].intersection(triple[2]) not in {frozenset(), triple[0], triple[2]}
            ):
                if triple[0].intersection(triple[1]) == triple[0].intersection(triple[2]):
                    return False
        return True
    
    def hasse_diagramm(self):
        G = nx.DiGraph()
        for pair in permutations(self.phyNet.clusters_by_node.values(),2):
            if pair[0].issubset(pair[1]) and pair[0] != pair [1]:
                skip = [pair[0], pair[1]]
                found_cluster_between = False
                for cluster in self.phyNet.clusters_by_node.values():
                    if cluster in skip:
                        continue
                    elif (cluster.issubset(pair[1]) and pair[0].issubset(cluster)):
                        found_cluster_between = True
                if not found_cluster_between:
                    G.add_edge(frozenset(pair[1]),frozenset(pair[0]))
        return G
        

    def _powerset(self, iterable, allow_empty_set:bool = True): 
        s = set(iterable)
        if allow_empty_set:
            return [set(subset) for subset in chain.from_iterable(combinations(s, r) for r in range(len(s)+1))]
        else: 
            return [set(subset) for subset in chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1))]





