"""Test prepReader Library

python3 test_prepReader_basic.py

In this cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.

"""

import sys
import unittest

sys.path.append('../')

try:
    import prepreader as pr
except ImportError:
    pr = None


class prepReaderLibrary(unittest.TestCase):
    def setUp(self):
        self.assertTrue(pr is not None, "Failed to import the prepReader module")

class TestprepReaderLibrary(prepReaderLibrary):

    def test_list_difference(self):

        """
        Fail if the provided list has elements of the first list which are in the
        second list.
        
        """
        list_1 = [1,2,3,4,5,6,7,8,9,10]
        list_2 = [8,9,10,11,12,13,14,15]
        self.assertEqual(pr.list_difference(list_1,list_2), [1,2,3,4,5,6,7])

    def test_uniquify(self):
    	
        """
        Fail if the provided list has dublicates.
        
        """

        list_1 = [0.4,-0.4,50,7,0,7]
        self.assertEqual(pr.uniquify(list_1), [0.4,-0.4,50,7,0])

    def test_getAtNum(self):

    	atoms = ["C","H","O","S"]
    	for element in atoms:
    		if element == "C":
    			# Fail the provided value is 6
    		    self.assertNotEqual(pr.getAtNum(element),1 )
    		elif element == "H":
    			# Fail the provided value is != 1
    		    self.assertEqual(pr.getAtNum(element),1 )
    		elif element == "O":
    			# Fail the provided value is = int (8) or = float(8)
    		    self.assertNotEqual(pr.getAtNum(element),8.0000001)
    		    # Fail the provided value is != int (8) or != float(8)
    		    self.assertEqual(pr.getAtNum(element),8.000000)



if __name__ == "__main__":
    unittest.main()