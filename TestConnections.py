# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import unittest
from src import rna2d

class  NewTestCase(unittest.TestCase):
    
    def test_connections_handling(self):
        assert rna2d.connection_dict_creator("....(((..[[...))).]]") == {0: None, 1: None, 2: None, 3: None, 4: 16, 5: 15, 6: 14, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: 6, 15: 5, 16: 4, 17: None, 18: None, 19: None}

    

if __name__ == '__main__':
    unittest.main()

