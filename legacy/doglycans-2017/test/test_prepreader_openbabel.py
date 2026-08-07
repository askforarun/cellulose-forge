"""Test prepReader-openbabel Function 

python3 test_prepReader_openbabel.py

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

try:
    import openbabel as ob
except ImportError:
    ob = None

try:
    import pybel
except ImportError:
    pybel = None


class prepReaderOpenbabel(unittest.TestCase):
    def setUp(self):
        self.assertTrue(pr is not None, "Failed to import the prepReader module")
        self.assertTrue(ob is not None, "Failed to import the openbabel module" )


class TestprepReaderOpenbabel(prepReaderOpenbabel):

    def test_getBonds(self):

        """
        Fail for different bond list 
        
        """
        pybmol = (pybel.readstring("smi", "CCCCCCCC"))
        mol = pybmol.OBMol       

        self.assertEqual(pr.getBonds(mol), [[1, 2], [2, 3], 
        	                                [3, 4], [4, 5], 
        	                                [5, 6], [6, 7], 
        	                                        [7, 8]])

    def test_getAngles(self):
    	
        """
        Fail for different number of angles 
        
        """
        pybmol = (pybel.readstring("smi", "CCCCCCCC"))
        mol = pybmol.OBMol       

        self.assertTrue(len(pr.getBonds(mol)),6)


    def test_getTorsion(self):
    	
        """
        Fail if not empty list
        
        """
        pybmol = (pybel.readstring("smi", "CC(=O)Cl"))
        mol = pybmol.OBMol       

        self.assertEqual(pr.getTorsions(mol), [] )


class PybelWrapper(prepReaderOpenbabel):
    def testDummy(self):
        self.assertTrue(pybel is not None, "Failed to import the Pybel module")


if __name__ == "__main__":
    unittest.main()