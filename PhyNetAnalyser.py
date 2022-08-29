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
        
        

    def is_phylogenetic(self):
        if self.phyNet.root is None: 
            return False
        for node in self.phyNet.digraph.nodes:
            if self.phyNet.digraph.out_degree(node) == 1 and self.phyNet.digraph.in_degree(node) <= 1:
                return False
        return True


    def is_level_k(self):
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
        return max_hybrids



    def is_tree_child(self):
        if self.is_tree_child is not None: 
            return self.is_tree_child

        tree_child = True
        self.is_tree_child = False
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
        return tree_child

    def is_pcc(self):
        # check if pcc property has already been calculated 
        if self.is_pcc is not None: 
            return self.is_pcc

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

        for node in self.phyNet.digraph.nodes:
            # iterate over hybrid nodes only
            if self.phyNet.digraph.in_degree(node) >= 2:
                predecessors = self.phyNet.digraph.predecessors(node)
                for predecessor in predecessors:
                    for child in self.phyNet.digraph.successors(predecessor):
                    # check for every predecessor of the hybrid node if there is a child that is not the hybrid node itself
                    # but there is also a path from the child of the predecessor to the hybrid node
                        if child is not node and nx.has_path(self.phyNet.digraph, child, node) :
                            self.is_shortcut_free = False
                            return False     
        self.is_shortcut_free = True                   
        return True 
        

    def is_normal(self,G):
        if self.is_shortcut_free() and self.is_tree_child():
            return True
        else:
            return False


    def is_semi_regular(self):
        if self.is_shortcut_free() and self.is_pcc():
            return True
        else:
            return False
    
    def is_regular(self):
        #semiregular kein knoten mit outdegree 1 
        return


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

        
    # lemma 3.30
    
    def is_closed(self):
        clustering_system = set([frozenset(cluster) for cluster in self.phyNet.clusters_by_node.values()])
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





