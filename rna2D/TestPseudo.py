# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


import unittest
from src import rna2d

class  NewTestCase(unittest.TestCase):
    
    def test_pseudo(self):
        assert rna2d.pseudo_handler("--<.(((..{[...)))>]}")=={0: None, 1: None, 2: 17, 3: None, 4: None, 5: None, 6: None, 7: None, 8: None, 9: 19, 10: 18, 11: None, 12: None, 13: None, 14: None, 15: None, 16: None, 17: 2, 18: 10, 19: 9}
	

if __name__ == '__main__':
    unittest.main()
	