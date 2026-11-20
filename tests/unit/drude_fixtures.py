"""Synthetic CHARMM Drude stream fixtures shared by the Drude unit tests."""

from __future__ import annotations

# A small invented residue and linkage that use every directive the parser
# supports.  It is not the official stream, which is not redistributed here.
MAIN = """\
* synthetic main stream
MASS -1 TC 12.01100 C
MASS -1 TO 15.99940 O
MASS -1 TE 15.99940 O
MASS -1 TH 1.00800 H
MASS -1 LPD 0.00000 H
MASS -1 DRUD 0.00000 H
read para card append
BONDS
TC TO 300.0 1.42
TO TH 500.0 0.96
TC TC 220.0 1.52
TC TH 300.0 1.10
TC TE 300.0 1.42
X DRUD 500.0 0.0
ANGLES
TO TC TC 50.0 110.0
TO TC TH 50.0 110.0
TC TC TH 50.0 110.0
TC TO TH 50.0 108.0
TC TE TC 50.0 115.0
TE TC TC 50.0 110.0
TE TC TH 50.0 110.0
DIHEDRALS
X TC TC X 0.2 3 0.0
X TC TO X 0.2 3 0.0
X TC TE X 0.2 3 0.0
NONBONDED nbxmod 5 atom vatom cdiel vdistance switch vswitch -
cutnb 16.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
TC 0.0 -0.05 2.0
END
"""

RESIDUE = """\
RESI BGLC 0.000 ! synthetic
GROUP
ATOM C1 TC 0.100 ALPHA -1.000 THOLE 1.300
ATOM H1 TH 0.100
ATOM O1 TO 0.000 ALPHA -0.800 THOLE 0.500
ATOM HO1 TH 0.300
ATOM LP1A LPD -0.250
ATOM LP1B LPD -0.250
GROUP
ATOM C4 TC 0.100 ALPHA -1.100 THOLE 0.600
ATOM H4 TH 0.100
ATOM O4 TO 0.000 ALPHA -0.900 THOLE 2.000
ATOM HO4 TH 0.300
ATOM LP4A LPD -0.250
ATOM LP4B LPD -0.250
BOND C1 O1 O1 HO1 C1 H1 C1 C4
BOND C4 O4 O4 HO4 C4 H4
BOND O1 LP1A O1 LP1B
BOND O4 LP4A O4 LP4B
LONEPAIR relative LP1A O1 C1 HO1 distance 0.35 angle 110.0 dihe 91.0
LONEPAIR relative LP1B O1 C1 HO1 distance 0.35 angle 110.0 dihe 269.0
ANISOTROPY O1 C1 LP1A LP1B A11 0.8108 A22 1.2162
LONEPAIR relative LP4A O4 C4 HO4 distance 0.35 angle 110.0 dihe 91.0
LONEPAIR relative LP4B O4 C4 HO4 distance 0.35 angle 110.0 dihe 269.0
ANISOTROPY O4 C4 LP4A LP4B A11 0.8108 A22 1.2162
IC O1 C4 *C1 H1 1.39 110.9 120.1 104.6 1.08
PATCH FIRST NONE LAST NONE
"""

PATCH = """\
PRES 14bb 0.000 ! synthetic
dele atom 1HO4
dele atom 1LP4A
dele atom 1LP4B
dele atom 2HO1
dele atom 2O1
dele atom 2DO1
dele atom 2LP1A
dele atom 2LP1B
ATOM 1C4 TC 0.050 ALPHA -1.050 THOLE 0.200
ATOM 1O4 TE 0.000 ALPHA -0.700 THOLE 1.000
ATOM 2C1 TC 0.150 ALPHA -1.020 THOLE 0.460
ATOM 1LPDA4 LPD -0.200
ATOM 1LPDB4 LPD -0.200
ATOM 1LPXD4 LPD 0.000 !dummy
BOND 1O4 2C1
BOND 1O4 1LPDA4 1O4 1LPDB4 1O4 1LPXD4
LONEPAIR bisector 1LPDA4 1O4 1C4 2C1 distance 0.35 angle 110.0 dihe 90.0
LONEPAIR bisector 1LPDB4 1O4 1C4 2C1 distance 0.35 angle 110.0 dihe 270.0
LONEPAIR bisector 1LPXD4 1O4 1C4 2C1 distance 0.10 angle 0.0 dihe 0.0
ANISOTROPY 1O4 1LPXD4 1LPDA4 1LPDB4 A11 0.8889 A22 1.2222
IC 1C4 1O4 2C1 2H1 1.45 121.0 -130.9 108.6 1.44
"""

TRAILER = "PRES OTHER 0.000\nend\n"


def carbohydrate(residue: str = RESIDUE, patch: str = PATCH) -> str:
    return "read rtf card append\n" + residue + "\n" + patch + "\n" + TRAILER
