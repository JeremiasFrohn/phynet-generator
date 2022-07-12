import networkx as nx


class PhyNetAnalyser: 
    def __init__(self, G:nx.DiGraph):
        #self.tree_child = None
        #self.shortcut_free = None
        self.digraph = G
        self.node_dict = self._create_node_dict()
        self.cluster_dict = self._create_cluster_dict()
        self.root = None
        

    def is_level_k(self, G:nx.DiGraph):
        max_hybrids = 0
        for biconnect_component in nx.biconnected_components(nx.Graph(G.to_undirected())):
            current_hybrid_counter = 0
            for node in biconnect_component:
                if self.digraph.in_degree(node) >= 2 and not set(G.predecessors(node)).isdisjoint(biconnect_component):
                    #ignore highest hybrid 
                    #if any(n in set(G.predecessors(node)) for n in biconnect_component)
                    current_hybrid_counter = current_hybrid_counter + 1
            if current_hybrid_counter > max_hybrids:
                max_hybrids = current_hybrid_counter
        return max_hybrids


    def is_tree_child(self, G:nx.DiGraph):
        tree_child = None
        for node in G.nodes:
            # iterate over every node that is not a leave 
            if G.out_degree(node) > 0:
                # set tree child false until at least one tree child could be found
                tree_child = False
                for child in G.successors(node):
                    if G.in_degree(child) == 1:
                        tree_child = True
                        break
                if not tree_child:
                    break
        return tree_child

    def is_pcc(self, G:nx.DiGraph):
        #pairwise iteration
        #itertools
        for node1 in G.nodes:
            for node2 in G.nodes: 
                if node2 == node1:
                    continue
                # check if nodes are comparable

                if node1 in self.node_dict[node2] or node2 in self.node_dict[node1]:
                    # if the nodes are comparable then check for subset condition
                    if self.cluster_dict[node1] <= self.cluster_dict[node2] or self.cluster_dict[node2] <= self.cluster_dict[node1]:
                        # if comparable and subset continue
                        continue
                    else:
                        # if comparable but not subset graph is no pcc
                        print(self.node_dict)
                        print("comparable nodes where subset condition fails")
                        print([node1,node2])
                        return False
                else: 
                    # if nodes are not comparable check for subset
                    if self.cluster_dict[node1] <= self.cluster_dict[node2] or self.cluster_dict[node2] <= self.cluster_dict[node1]:
                        # if subset but not comparable graph is no pcc else continue 
                        print(self.node_dict)
                        print("non comparable nodes where subset condition is true")
                        print([node1,node2])
                        return False
        return True 
    
    def is_shortcut_free(self, G:nx.DiGraph):
        root = self._root(G)
        for node in G.nodes:
            if G.in_degree(node) >= 2:
                if len([path for path in nx.all_shortest_paths(G, root, node)]) < 2:
                    return False
        return True 
        

    def is_normal(self,G):
        if self.is_shortcut_free(G) and self.is_tree_child(G):
            return True
        else:
            return False

    
    def _root(self,G):
        for node in G.nodes:
            if G.in_degree(node) == 0:
                return node
        return None

    def semi_regular():
        return


    def _reachable_nodes(self, n):
        # list of nodes that need to be checked for succcessors
        nodes_to_check = [successor for successor in self.digraph.successors(n)]
        nodes_to_check_copy = list(nodes_to_check)
        # list to return
        reachable_nodes = set(nodes_to_check)
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


    def _create_node_dict(self):
        node_dict = dict()
        for node in self.digraph.nodes:
            node_dict[node] = self._reachable_nodes(node)
        return node_dict
    

    def _create_cluster_dict(self):
        cluster_dict = dict()
        for node, reachables in self.node_dict.items():
            reachable_leaves = set()
            for reachable in reachables: 
                if self.digraph.out_degree(reachable) == 0:
                    reachable_leaves.add(reachable)
            cluster_dict[node] = reachable_leaves
        return cluster_dict
    

