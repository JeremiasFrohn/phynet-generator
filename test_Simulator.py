import unittest
import Simulator as sim
import networkx as nx

class Test_Simulator(unittest.TestCase):

    def setUp(self):
        self.sim1 = sim.Simulator()
        self.G = nx.DiGraph()
        self.G.add_edges_from([(1, 2), (1, 3), (2, 4), (2, 5), (4, 6), (4, 7), (6, 8), (6, 7), (7, 9), (8, 10), (8, 11)])
        self.sim1.digraph = self.G


    def test_reachable_nodes(self):

        self.assertEqual(self.sim1.reachable_nodes(1), set([2,3,4,5,6,7,8,9,10,11]))
        
        with self.assertRaises(nx.NetworkXError):
            self.sim1.reachable_nodes(12)

        self.assertEqual(self.sim1.reachable_nodes(11), set())

if __name__ == '__main__':
    unittest.main()
