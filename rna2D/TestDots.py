# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


import unittest
from src import rna2d

class  NewTestCase(unittest.TestCase):
    
    def test_dots(self):
        assert rna2d.dots_handler("(((((((.........(((........))))))))))......((((..........))))...........",{0: 36, 1: 35, 2: 34, 3: 33, 4: 32, 5: 31, 6: 30, 7: None, 8: None, 9: None, 10: None, 11: None, 12: None, 13: None, 14: None, 15: None, 16: 29, 17: 28, 18: 27, 19: None, 20: None, 21: None, 22: None, 23: None, 24: None, 25: None, 26: None, 27: 18, 28: 17, 29: 16, 30: 6, 31: 5, 32: 4, 33: 3, 34: 2, 35: 1, 36: 0, 37: None, 38: None, 39: None, 40: None, 41: None, 42: None, 43: 60, 44: 59, 45: 58, 46: 57, 47: None, 48: None, 49: None, 50: None, 51: None, 52: None, 53: None, 54: None, 55: None, 56: None, 57: 46, 58: 45, 59: 44, 60: 43, 61: None, 62: None, 63: None, 64: None, 65: None, 66: None, 67: None, 68: None, 69: None, 70: None, 71: None}) == [[(61, 71)], [(19, 26), (47, 56)], [], [(7, 15)], [(37, 42)]]

    

if __name__ == '__main__':
    unittest.main()