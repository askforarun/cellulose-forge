"""Test prepReader-openbabel Function 

python3 test_seqparse.py

In this cases, the test file is run directly from the source folder,
and so you can quickly develop the tests and try them out.

"""

import sys
import unittest

sys.path.append('../')

try:
    import seqparse
except ImportError:
    seqparse = None


class prepReaderInputReader(unittest.TestCase):
    def setUp(self):
        self.assertTrue(seqparse is not None, "Failed to import the seqparse module")


class TestprepReaderOpenbabel(prepReaderInputReader):

    def test_seqFileReader(self):


        """
        Fail if the provided lists have different elements from the output.
        
        """
        s1 = "HM=P.150/HD22:-(ND2,C1)-{4YB-(O4,C1)-4YB<C3,C4>-(O4,C1,<Con,C5>)}3" \
        "-VMB-(O6,C1,[H5 C5 C6 O6 0.0 C5 C6 O6 C1 180.0 C6 O6 C1 H1 60.0])-0MA"

        self.assertEqual(seqparse.parser.parse(s1)[0],  ("HM"))
        self.assertEqual(seqparse.parser.parse(s1)[1],  ('P.150', 'HD22'))
        self.assertEqual(seqparse.parser.parse(s1)[2],  (('ND2', 'C1'), None ))
        self.assertEqual(seqparse.parser.parse(s1)[3],  ("4YB"))
        self.assertEqual(seqparse.parser.parse(s1)[4],  (('O4', 'C1'), None))
        self.assertEqual(seqparse.parser.parse(s1)[5],  ("4YB"))
        self.assertEqual(seqparse.parser.parse(s1)[6],  (('O4', 'C1'), 'Con', 'C3', 'C4'))
        self.assertEqual(seqparse.parser.parse(s1)[7],  ("4YB"))
        self.assertEqual(seqparse.parser.parse(s1)[8],  (('O4', 'C1'), None))
        self.assertEqual(seqparse.parser.parse(s1)[9],  ("4YB"))
        self.assertEqual(seqparse.parser.parse(s1)[10],(('O4', 'C1'), 'Con', 'C3', 'C4'))
        self.assertEqual(seqparse.parser.parse(s1)[11],("4YB"))
        self.assertEqual(seqparse.parser.parse(s1)[12],(('O4', 'C1'), None))
        self.assertEqual(seqparse.parser.parse(s1)[13],("4YB"))
        self.assertEqual(seqparse.parser.parse(s1)[14],(('O4', 'C1'), 'Con', 'C3', 'C4'))
        self.assertEqual(seqparse.parser.parse(s1)[15],("VMB"))
        self.assertEqual(seqparse.parser.parse(s1)[16],(('O6', 'C1'), [('H5', 'C5', 'C6', 'O6', 0.0), ('C5', 'C6', 'O6', 'C1', 180.0), ('C6', 'O6', 'C1', 'H1', 60.0)]))
        self.assertEqual(seqparse.parser.parse(s1)[17],("0MA"))


if __name__ == "__main__":
    unittest.main()