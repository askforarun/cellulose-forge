#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Version 1
Author: See the AUTHORS file for a list of contributors to doGlycans.py
License GPLv3

Reads Glycam PREP database and builds saccharides specified by sequence file.
Writes structure to coordinate file and Gromacs topology to ITP file.

"""

import re
import sys
import random
import os.path
import openbabel as ob
from collections import OrderedDict
from optparse import OptionParser
from seqparse import parser
from math import sqrt


isymdu = "DU"


def getFNs(option, opt_str, value, parser):
    """Extract number of parameters from command line.

    """
    reX = re.compile(r'^[^-]+\S+$')
    value = []
    rargs = parser.rargs
    while rargs:
        arg = rargs[0]
        if reX.match(arg) is None:
            break
        else:
            value.append(arg)
            del rargs[0]
    setattr(parser.values, option.dest, value)


def optP():
    """Parse command line options.

    """
    usage = "[python3] %prog [-Wmdt] -f $DGPATH/dat/glycam06h.itp" \
            "-p $DGPATH/dat/GLYCAM_06h.prep -s sequence_file.seq" \
            "-c POL.pdb -o POL.itp --opls --amber --charmm"

    description = "Reads Glycam PREP database and builds saccharides" \
                  "specified by sequence file. Writes structure" \
                  "to coordinate file and Gromacs topology to ITP file."
    version = "\n%prog Version 0.02\n\n" \
              "Requires Python 3, Open Babel Python bindings and PLY."

    optParser = OptionParser(usage=usage,
                             version=version,
                             description=description)

    optParser.set_defaults(ffFN=None, prepFNs=['Glycam_06.prep'], seqFN=None,
                           itpFN=None, molFN=None, opls=True, charmm=False,
                           amber=False, gen_mono=False, gen_di=False, 
                           gen_tri=False, no_warn=False)

    optParser.add_option('-f', type='str',
                         dest='ffFN',
                         help="Glycam parameter ITP file (in)"
                         " [default: %default]")
    optParser.add_option('-p', action='callback', callback=getFNs,
                         dest='prepFNs',
                         help="Prep file(s) (in)"
                         " [default: %default]")
    optParser.add_option('-s', type='str',
                         dest='seqFN',
                         help="Sequence file (in)"
                         " [default: %default]")
    optParser.add_option('-m', action='store_true', default=False,
                         dest='gen_mono',
                         help="Generate PDBs and ITPs for monomers"
                         " [default: %default]")
    optParser.add_option('-d', action='store_true',
                         dest='gen_di',
                         help="Generate PDBs and ITPs for dimers"
                         " [default: %default]")
    optParser.add_option('-t', action='store_true',
                         dest='gen_tri',
                         help="Generate PDBs and ITPs for trimers"
                         " [default: %default]")
    optParser.add_option('-W', action='store_true', default=False,
                         dest='no_warn',
                         help="Suppress warnings"
                         " [default: %default]")
    optParser.add_option('-c', type='str',
                         dest='molFN',
                         help="Coordinate file (out)"
                         " [default: %default]")
    optParser.add_option('-o', type='str',
                         dest='itpFN',
                         help="Gromacs topology file (out)"
                         " [default: %default]")
    optParser.add_option('--opls', action='store_true',
                         dest='opls',
                         help="Add suffix '_S' for opls atom types"
                         " [default: %default]")
    optParser.add_option('--charmm', action='store_true',
                         dest='pathFN',
                         help="Print PDB for charmm-gui submission."
                         " [default: %default]")
    optParser.add_option('--amber', action='store_true',
                         dest='amber',
                         help="Remove suffix '_S' for amber atom types"
                         " [default: %default]")

    return optParser.parse_args()


class SeqRes(object):
    """Struct for storing sequence residue.

    Parameters can be obtained using attribute access syntax.
    """
    def __init__(self, **kwds):
        """Initialize default propreties and set values from arguments.

        """
        defaults = {k: None for k in ('name',
                                      'trunk_id',
                                      'branch_id',
                                      'bond',
                                      'dihedrals',
                                      'replace')}
        self.__dict__.update(defaults)
        self.__dict__.update(kwds)

    def __str__(self):
        """Make string representation of the object for printing.

        """
        s = ""
        defaults = ['name',
                    'trunk_id',
                    'branch_id',
                    'bond',
                    'dihedrals',
                    'replace']
        for k in defaults:
            s += "{0:12s} {1}\n".format(k+':', self.__dict__[k])
        for k, v in self.__dict__.items():
            if k not in defaults:
                s += "{0:12s} {1}\n".format(k+':', v)

        return s


def getAtNum(atype):
    """Deduce atom number from the first letter of atomtype name.

    """
    if atype[0] == 'C':
        return 6
    elif atype[0] == 'H':
        return 1
    elif atype[0] == 'O':
        return 8
    elif atype[0] == 'N':
        return 7
    elif atype[0] == 'S':
        return 16
    elif atype[0] == 'F':
        if len(atype) == 1:
            return 9
        elif atype[1] in 'Ee':
            return 26
    elif atype == isymdu:
        return 0
    else:
        print("Error: Unrecognized element for atomtype", atype)


def getMolAtomByName(mol, aname):
    """Returns atom from molecule with specified name (IGRAPH).

    """
    for i in ob.OBMolAtomIter(mol):
        if i.IGRAPH == aname:
            return i


def getResAtomByName(res, aname):
    """Returns atom from residue with specified name.

    """
    for i in ob.OBResidueAtomIter(res):
        if res.GetAtomID(i) == aname:
            return i

    s = "Error: Could not find atom {} in residue {}"
    print(s.format(aname, res.GetData('BRANCH_ID').GetValue()))
    for i in ob.OBResidueAtomIter(res):
        print(res.GetAtomID(i))
    raise KeyError


def getResByBranchID(mol, branch_id):
    """Return residue with given branch id.

    """
    for i in ob.OBResidueIter(mol):
        if i.GetData("BRANCH_ID").GetValue() == branch_id:
            return i

    s = "Error: Could not find residue with BRANCH_ID={}"
    print(s.format(branch_id))
    raise KeyError


def getAtomNbrByName(atom, aname):
    """Return atom with specified name bonded to given atom.

    """
    for i in ob.OBAtomAtomIter(atom):
        if i.HasData("IGRAPH"):
            if i.IGRAPH == aname:
                return i
        else:
            if i.GetResidue().GetAtomID(i).strip() == aname:
                return i

    s = "Error: Could not find neighbor named {} for atom {}"
    if atom.HasData("IGRAPH"):
        print(s.format(aname, atom.IGRAPH))
    else:
        print(s.format(aname, atom.GetType()))
    for i in ob.OBAtomAtomIter(atom):
        if i.HasData("IGRAPH"):
            print(i.IGRAPH)
        else:
            print(i.GetResidue().GetAtomID(i))
    raise KeyError


def addBond(mol, begin, end, bond=ob.OBBond()):
    """Add to molecule a bond with given start and end atoms.

    """
    bond.SetBegin(begin)
    bond.SetEnd(end)
    return mol.AddBond(bond)


def addPairData(mol, name, value, d=ob.OBPairData()):
    """Add OBPairData to the unit (atom, molecule or residue).

    """
    d.SetAttribute(name)
    d.SetValue(value)
    mol.CloneData(d)


def mkProperty(cls, name):
    """Make property and put it in given class as attribute with given name.

    It allows to add and get named PairData using attribute reference syntax.
    """
    def getProperty(self):
        return self.GetData(name).GetValue()

    def setProperty(self, value):
        if self.HasData(name):
            self.DeleteData(name)
        addPairData(self, name, value)

    def delProperty(self, name):
        self.DeleteData(name)

    setattr(cls, name, property(getProperty, setProperty, delProperty, name))

# Set OBAtom properties for information stored in PairData
for n in ['I', 'IGRAPH', 'ISYMBL', 'ITREE', 'NA', 'NB', 'NC',
          'R', 'THETA', 'PHI', 'CHRG']:
    mkProperty(ob.OBAtom, n)


def getImproper(mol, impStr):
    """Finds improper torsion atom indices.

    Takes a OBMol instance and dash-seperated string of atom names as arguments
    and returns a list of corresponding atom indices.
    """
    impNames = impStr.split('-')
    imp = [getAtomByName(mol, i).GetIdx() for i in impNames]

    return imp


def readPrepFile(f, no_warn=False):
    """The main function for reading prep files.

    As an argument takes file handler and returns a dictionary with the residue
    names as keys and OBMol objects as values.
    Note: All dummy atoms in molecules are removed.
    """
    # {'Residue name': [OBMol, ...extra information]}
    residues = {}

    line = f.readline().split()

    try:
        if line[2] != '2':
            print("Error: Unknown database type", line[2])
            sys.exit(2)
    except IndexError:
        # Missing PREP file header
        pass

    f.readline()

    while line != '':
        title = f.readline()
        if title.strip() == 'STOP':
            # End of the prep file
            break
        f.readline()
        line = f.readline().split()
        resname = line[0]
        try:
            # Flag for the type of coordinates to be saved for the LINK module.
            # We are not using it.
            intx = line[1]
            if line[2] != '0':
                print("Error: Supported are only formatted files!")
                sys.exit(2)
        except IndexError:
            intx = "INT"

        line = f.readline().split()
        ifixc = line[0]
        if ifixc != 'CORRECT':
            print("Error: Not implemented IFIXC=" + ifixc)
            sys.exit(5)
        iomit = line[1]
        if iomit != 'OMIT':
            print("Error: Not implemented IOMIT=" + iomit)
            sys.exit(5)
        global isymdu
        isymdu = line[2]
        if isymdu != 'DU':
            print("Warning: It is preferable to use 'DU' for dummy atoms!")
        # Which dummy atoms to delete.
        # Open Babel removes all dummy atoms in InternalToCartesian() and
        # there is no option to disable it.
        # Ignoring this setting.
        ipos = line[3]

        cut = float(f.readline())
        if cut != 0.0:
            # TODO: Calculate loop closing bonds
            #print("Warning: Non-zero cutoff for loop closing bonds ignored!")
            pass

        # Create an empty molecule
        mol = ob.OBMol()
        mol.SetTitle(title.strip('\n'))
        mol.SetAutomaticPartialCharge(False)

        # Create residue in the molecule
        res = mol.NewResidue()
        res.SetName(resname)
        addPairData(res, 'NOTE', title.strip(' \n'))

        # Read atoms
        line = f.readline()
        while line.strip() != '':
            line = line.split()

            # Create and populate the OBAtom object
            atom = mol.NewAtom()
            atom.SetAtomicNum(getAtNum(line[2]))
            atom.SetPartialCharge(float(line[10]))

            #atom.I      = line[0]
            atom.IGRAPH = line[1]
            atom.ISYMBL = line[2]
            #atom.ITREE  = line[3]
            atom.NA     = line[4]
            atom.NB     = line[5]
            atom.NC     = line[6]
            atom.R      = line[7]
            atom.THETA  = line[8]
            atom.PHI    = line[9]
            #atom.CHRG   = line[10]

            if line[2] != isymdu:
                res.AddAtom(atom)
                res.SetAtomID(atom, line[1])

            # Read the next line
            line = f.readline()

        # Make Z-matrix
        internalCoords = []
        # The first element in internalCoords has to be NULL for old Open Babel
        if mol.NumAtoms() + 1 == len(mol.GetInternalCoord()):
            internalCoords.append(None)
        for obaa in ob.OBMolAtomIter(mol):
            na = int(obaa.NA)
            nb = int(obaa.NB)
            nc = int(obaa.NC)
            if na == 0:
                # The first atom
                if obaa.ISYMBL != isymdu:
                    print("Error: The first three atoms should be", end=' ')
                    print(isymdu, "in", resname)
                    sys.exit(2)
                ic = ob.OBInternalCoord()
            elif nb < 1:
                # The second atom
                if obaa.ISYMBL != isymdu:
                    print("Error: The first three atoms should be", end=' ')
                    print(isymdu, "in", resname)
                    sys.exit(2)
                r = float(obaa.R)
                A = mol.GetAtom(na)
                addBond(mol, A, obaa)
                ic = ob.OBInternalCoord(a=A, dst=r)
            elif nc < 1:
                # The third atom
                if obaa.ISYMBL != isymdu:
                    print("Error: The first three atoms should be", end=' ')
                    print(isymdu, "in", resname)
                    sys.exit(2)
                r = float(obaa.R)
                theta = float(obaa.THETA)
                if theta == 0.0:
                    if not no_warn:
                        print("Warning: THETA=0.0 for the third atom", end=' ')
                        print("of ", resname, ", setting to 90.0", sep='')
                    theta = 90.0
                A = mol.GetAtom(na)
                B = mol.GetAtom(nb)
                addBond(mol, A, obaa)
                ic = ob.OBInternalCoord(a=A, b=B, dst=r, ang=theta)
            else:
                # The rest
                if obaa.ISYMBL == isymdu:
                    print("Warning: All dummy atoms will be removed", end=' ')
                    print("from residue", resname)
                r = float(obaa.R)
                theta = float(obaa.THETA)
                phi = float(obaa.PHI)
                A = mol.GetAtom(na)
                B = mol.GetAtom(nb)
                C = mol.GetAtom(nc)
                addBond(mol, A, obaa)
                ic = ob.OBInternalCoord(a=A, b=B, c=C,
                                        dst=r, ang=theta, tor=phi)
            internalCoords.append(ic)

        # Set coordinates for the molecule
        # For this to work, have to apply a patch for openbabel-python (r4400)
        ob.InternalToCartesian(ob.vectorpOBInternalCoord(internalCoords), mol)

        # Dummy atoms are removed by ob.InternalToCartesian(), but a patch has
        # to be applied to OB that fixes that otherwise only every other dummy
        # atom is removed.

        # This removes the one dummy atom left if OB wasn't patched.
        for a in ob.OBMolAtomIter(mol):
            if a.ISYMBL == isymdu:
                mol.DeleteAtom(a)
            break

        iopr = f.readline()
        while iopr.strip() != 'DONE':
            iopr = iopr.strip()
            if iopr == 'CHARGE':
                # Read adittional charges
                print("Error: Not implemented IOPR=" + iopr)
                sys.exit(5)
            elif iopr == 'LOOP':
                # Read explicit loop closing bonds
                lbonds = []
                line = f.readline()
                while line.strip() != '':
                    line = line.split()
                    lbonds.append(line[:3])
                    line = f.readline()
                for b in lbonds:
                    addBond(mol,
                            getMolAtomByName(mol, b[0]),
                            getMolAtomByName(mol, b[1]))
                iopr = f.readline()
                continue
            elif iopr == 'IMPROPER':
                # Read improper torsion angles
                impropers = []
                line = f.readline()
                while line.strip() != '':
                    line = line.split()
                    impropers.append(line[:5])
                    line = f.readline()
                impStr = ""
                for i in impropers:
                    a = i[0]
                    b = i[1]
                    c = i[2]
                    d = i[3]
                    impStr += "{0:s}-{1:s}-{2:s}-{3:s} ".format(a, b, c, d)
                # Skip impropers on anomeric carbon
                if impropers != [] and b != 'H1':
                    addPairData(mol, 'IMPROPER', impStr.strip())
                iopr = f.readline()
                continue
            elif iopr == '':
                # Got more than one blank line
                iopr = f.readline()
            else:
                print("Error: Unsupported IOPR=" + iopr)
                sys.exit(2)

        residues[resname] = mol

    return residues


def unpackConnection(c):
    """Unpack connection data coming from parser.

    """
    # connection = ((bond), [angles], ['branches'])
    if c[0] is not None:
        con = [c[0][0], c[0][1]]
    else:
        # Terminal label is packed as connection without bond and angles
        con = [None, None]
    if c[1] is not None:
        con.append([angle for angles in c[1]
                              for angle in angles])
    return con


def getConnection(sres, branch, labels):
    """Add connection data to SeqRes object.

    """
    if isinstance(branch[0][1], str):
        # It is a label for angles, expand it
        sres.dihedrals = labels[branch[0][1]][2]
        sres.bond      = labels[branch[0][1]][:2]
    else:
        # It is a connection with or without angles
        c = unpackConnection(branch[0])
        sres.bond = c[:2]
        if len(c) > 2:
            sres.dihedrals = c[2]


def getAttachments(branch, branches, attachments, bonds, labels, loops = {}):
    """Recursive generator for residues of a given branch and its branches.

    Returns a residue name and a list of trunk id, branch id, the two atoms
    forming the bond between residues and optional torsion specified by four
    atom names and a value.
    """
    batts = {k: v for k, v in attachments.items()
                      if k.split('.')[0] == branch[0].branch_id.split('.')[0]}
    for attachment in batts.keys():
        for bbranch in batts[attachment]:
            # Get other end of the bond
            if bbranch in bonds:
                bonds[bbranch][1] = attachment
                continue
            # Update connection
            try:
                branches[bbranch][0].trunk_id = attachment
                for i in branches[bbranch]:
                    yield i
                # Recurse
                for i in getAttachments(branches[bbranch],
                                        branches, attachments,
                                        bonds, labels, loops):
                    yield i
            except KeyError:
                # This is not a branch, must be a loop
                L, a = bbranch.split('.')
                if L not in loops:
                    # Not seen, save it for later
                    loops[L] = [a, attachment]
                    continue
                else:
                    # Make a bond between the two
                    c = SeqRes(trunk_id=loops[L][1],
                               branch_id=attachment,
                               bond=[loops[L][0], a])
                    if L in labels:
                        # Add angles from label to connection
                        # TODO: Ensure correct atom order
                        c.dihedrals = labels[L][2]
                    yield c


def getSeqRes(sequence):
    """Generator to get next residue name in the sequence.

    Returns a residue name and a list of trunk id, branch id, the two atoms
    forming the bond between residues and optional torsion specified by four
    atom names and a value.
    """
    # 1. Construct dictionary of labels
    # 2. Expand labels for angles
    # 3. Resolve attachment points
    # 4. Flatten the sequence
    s = "{}.{}"
    reXe = re.compile(r'^\s*$')
    reXc = re.compile(r'^\s*#')
    bonds = OrderedDict()
    branches = OrderedDict()
    attachments = OrderedDict()
    replaceableAtom = None
    labels = {}
    for i in sequence:
        # Skip empty lines and comments
        if (reXe.match(i) is not None
            or reXc.match(i) is not None):
            continue
        # branch = [name, attach, connection, residue, ...]
        # attach = None | (parent, None | replace)
        # connection = None | ((bond), None | [angles]{, 'branch'})
        # connection_label = [name, None, connection]
        branch = parser.parse(i)
        if len(branch) == 3:
            if branch[1] is None:
                # It is label for angles
                labels[branch[0]] = unpackConnection(branch[2])
            else:
                # It is an explicit bond
                bonds[branch[0]] = [branch[1][0], None]
                bonds[branch[0]].extend(unpackConnection(branch[2]))
        else:
            # It is branch
            branches[branch[0]] = branch
            # Collect attachment points
            count = 1
            for res, con in zip(branch[3:-1:2], branch[4:-1:2]):
                if len(con) > 2:
                    attachments[s.format(branch[0], count)] = con[2:]
                count += 1
            if branch[-1][0] is None:
                # Add terminal label
                attachments[s.format(branch[0], count)] = branch[-1][2:]

    for branch in branches.values():
        count = 1
        branch_name = branch[0]

        # Initial residue needs special treatment in case it is attached
        # to another sequence
        branch_id = s.format(branch_name, count)
        if branch[1] is not None:
            trunk_id = branch[1][0]
            # Check for replacable atom
            if branch[1][1] is not None:
                replaceableAtom = branch[1][1]
        else:
            trunk_id = s.format(branch_name, count-1)
        del branch[:2]

        # List of SeqRes objects
        tseq = []
        tseq.append(SeqRes(name=branch[1],
                           trunk_id=trunk_id,
                           branch_id=branch_id))

        # Check connection
        if branch[0] is not None:
            getConnection(tseq[-1], branch, labels)

        # Save replacable atom
        if replaceableAtom is not None:
            tseq[-1].replace = replaceableAtom
            replaceableAtom = None

        del branch[:2]
        count += 1

        # Add remaining residues to the list
        while branch:
            # Check for terminal tag
            if len(branch) == 1:
                # There is only a tag there
                del branch[0]
                break
            trunk_id = s.format(branch_name, count-1)
            branch_id = s.format(branch_name, count)

            tseq.append(SeqRes(name=branch[1],
                               trunk_id=trunk_id,
                               branch_id=branch_id))
            getConnection(tseq[-1], branch, labels)

            del branch[:2]
            count += 1

        branches[branch_name] = tseq

    # Return all trunks and their branches and loops
    for branch in branches.values():
        if (branch[0].trunk_id.split('.')[0] == 'P'
            or branch[0].trunk_id.split('.')[1] == '0'):
            for i in branch:
                yield i
            # Get all attachments for the current branch recursively
            for i in getAttachments(branch,
                                    branches, attachments,
                                    bonds, labels):
                yield i

    # Add explicit bonds
    for b in bonds.keys():
        yield SeqRes(bond=bonds[b])


def setTorsions(mol, pres, nres):
    """Set torsions to given values when connecting residues.

    """
    if nres.dihedrals is not None:
        for n in range(len(nres.dihedrals) // 5):
            a1 = getResAtomByName(pres, nres.dihedrals[0+n*5])
            a2 = getAtomNbrByName(a1,   nres.dihedrals[1+n*5])
            a3 = getAtomNbrByName(a2,   nres.dihedrals[2+n*5])
            a4 = getAtomNbrByName(a3,   nres.dihedrals[3+n*5])
            angle = float(nres.dihedrals[4+n*5]) * ob.DEG_TO_RAD
            mol.SetTorsion(a1, a2, a3, a4, angle)


def buildSeq(residues, sequence):
    """Build a molecule with given sequence.

    """
    molname = sequence[0].split('=')[0]
    builder = ob.OBBuilder()
    # Use ring conformations from existing 3D coordinates
    builder.SetKeepRings()

    print("Seq:")
    for i in sequence:
        print(i)

    # Get the first residue
    lres = getSeqRes(sequence)
    nres = next(lres)
    mol = ob.OBMol(residues[nres.name])
    mol.SetAutomaticPartialCharge(False)
    addPairData(mol, 'PartialCharges', "ESP")

    res = mol.GetResidue(mol.NumResidues()-1)
    res.SetNum(res.GetIdx()+1)
    res.CloneData(residues[nres.name].GetResidue(0).GetData('NOTE'))
    res.CloneData(residues[nres.name].GetData('IMPROPER'))
    addPairData(res, 'TRUNK_ID',  nres.trunk_id)
    addPairData(res, 'BRANCH_ID', nres.branch_id)

    for nres in lres:
        if nres.name is not None:
            # Add next residue
            mol += residues[nres.name]

            # Have to manually copy residue data, += operator doesn't do it
            res = mol.GetResidue(mol.NumResidues()-1)
            res.SetNum(res.GetIdx()+1)
            res.CloneData(residues[nres.name].GetResidue(0).GetData('NOTE'))
            res.CloneData(residues[nres.name].GetData('IMPROPER'))
            addPairData(res, 'TRUNK_ID',  nres.trunk_id)
            addPairData(res, 'BRANCH_ID', nres.branch_id)

            # Get parent residue from the right chain
            p = getResByBranchID(mol, nres.trunk_id)

            # Get indices for atoms in connection
            a = getResAtomByName(p,   nres.bond[0]).GetIdx()
            b = getResAtomByName(res, nres.bond[1]).GetIdx()

            # Join the residues
            if builder.Connect(mol, a, b):
                # Set torsion(s)
                setTorsions(mol, p, nres)
            else:
                print("Error: Something went wrong while joining fragments!")
                sys.exit(3)
        else:
            # Don't add new residue, only make given connection
            r = getResByBranchID(mol, nres.branch_id)

            # Get parent residue
            p = getResByBranchID(mol, nres.trunk_id)

            # Get atoms
            a = getResAtomByName(p, nres.bond[0])
            b = getResAtomByName(r, nres.bond[1])

            # Make the bond
            addBond(mol, a, b)

            # Set torsion(s)
            setTorsions(mol, p, nres)

    mol.SetTitle(molname)

    return mol


def uniquify(seq):
    """Take a list and return list without duplicates.

    """
    seen = set()
    return [i for i in seq if str(i) not in seen and not seen.add(str(i))]


def list_difference(s1, s2):
    """Take set difference of two lists.

    Return a list of only those elements of the first list which are not in the
    second list.
    """
    s2 = set([str(i) for i in s2])
    return [j for j in s1 if str(j) not in s2]


def getResCharge(residue):
    """Calculate summary charge of given residue.

    """
    charge = 0
    for i in ob.OBResidueAtomIter(residue):
        charge += i.GetPartialCharge()

    return charge


def get14NBParams(ffFN):
    """Get 1-4 non-bonded intraction parameters and scaling factors.

    """
    # Regular expressions used
    reXa = re.compile(r'^[^;]*\[\s+atomtypes\s+\]')
    reXt = re.compile(r'^[^;]*\[\s+dihedraltypes\s+\]')
    reXc = re.compile(r'^\s*[;]+')
    reXe = re.compile(r'^\s*$')
    reXs = re.compile(r'^[^;]*\[\s+[a-z]*\s+\]')

    # Read Glycam forcfield file and extract improper torsions
    sc14 = OrderedDict()
    vdw = OrderedDict()
    with open(ffFN, 'r') as f:
        # Find [ atomtypes ]
        line = f.readline()
        while reXa.match(line) is None:
            line = f.readline()
        line = f.readline()

        # Read atomtypes and save VdW parameters
        while line != '' and reXs.match(line) is None:
            if (reXc.match(line) is not None
                or reXe.match(line) is not None):
                line = f.readline()
                continue
            a, n, m, ch, ptype, sigma, epsilon = line.split()[:7]
            vdw[a] = (float(sigma), float(epsilon))
            line = f.readline()

        # Find [ dihedraltypes ]
        line = f.readline()
        while reXt.match(line) is None:
            line = f.readline()

        # Read scaling factors
        line = f.readline()
        while line != '' and reXs.match(line) is None:
            if (reXc.match(line) is not None
                or reXe.match(line) is not None):
                line = f.readline()
                continue
            try:
                a, b, c, d, fn, p, kd, pn, com, scee, scnb = line.split()[:11]
            except ValueError:
                line = f.readline()
                continue
            if scee[:5] != 'SCEE=' or scnb[:5] != 'SCNB=':
                line = f.readline()
                continue
            sc14[(a,b,c,d)] = (float(scee[5:]), float(scnb[5:]))
            line = f.readline()

    return sc14, vdw


def getImproperTypes(ffFN):
    """Get improper torsions from forcefield ITP file.

    """
    # Regular expressions used
    reXt = re.compile(r'^[^;]*\[\s+dihedraltypes\s+\]')
    reXc = re.compile(r'^\s*[;]+')
    reXe = re.compile(r'^\s*$')
    reXs = re.compile(r'^[^;]*\[\s+[a-z]*\s+\]')

    # Read Glycam forcfield file and extract improper torsions
    imps = OrderedDict()
    with open(ffFN, 'r') as f:
        # Find second [ dihedraltypes ]
        line = f.readline()
        while reXt.match(line) is None:
            line = f.readline()
        line = f.readline()
        while reXt.match(line) is None:
            line = f.readline()
        # Read atomtypes
        line = f.readline()
        while line != '' and reXs.match(line) is None:
            if (reXc.match(line) is not None
                or reXe.match(line) is not None):
                line = f.readline()
                continue
            a, b, c, d, func, phase, kd, pn = line.split()[:8]
            imps[(a,b,c,d)] = (int(func), float(phase), float(kd), int(pn))
            line = f.readline()

    return imps


def getOrderedImp(atoms, imp):
    """Reorder atoms and check if bonding corresponds to impropers.

    a-b-c-d
    a and b are not bonded
    a, b, d are bonded to c
    Reorder atoms as in imp.
    """
    # If there are no wildcards in imp then use explicit atom order
    if 'X' not in imp:
        # Check if matched the correct one by comparing atomtype counts
        ats = [a.ISYMBL for a in atoms]
        if {t: ats.count(t) for t in ats} != {t: imp.count(t) for t in imp}:
            return None

        # Reorder atoms as in the improper
        iatoms = []
        for at in imp:
            for atom in atoms:
                if atom.ISYMBL == at:
                    if atom not in iatoms:
                        iatoms.append(atom)
                        break
        atoms = iatoms

        # Temporarily save central atom and remove it from the list
        c = atoms[2]
        atoms.remove(c)

        # Sort atoms with identical types by their indices
        pats = list(imp)
        pats.remove(imp[2])
        for i in pats:
            if pats.count(i) == 2:
                # Find which are they and swap if necessary
                ats = [(k, a) for k, a in enumerate(atoms)
                                  if a.ISYMBL == i]
                va = [a.GetIdx() for k, a in ats]
                ka = [k for k, a in ats]
                if va[0] > va[1]:
                    atoms[ka[0]], atoms[ka[1]] = atoms[ka[1]], atoms[ka[0]]
                break
            elif pats.count(i) == 3:
                atoms.sort(key=lambda a: a.GetIdx())
                break
    else:
        # Temporarily save central atom and remove it from the list
        for atom in atoms:
            if atom.ISYMBL == imp[2]:
                    c = atom
                    break
        atoms.remove(c)

        # Sort explicit atoms
        eatoms = [a for a in atoms
                        if a.ISYMBL in imp]
        eatoms.sort(key=lambda a: (a.ISYMBL, a.GetIdx()))

        # Atoms matched by wildcard are in front and sorted by index
        watoms = [a for a in atoms
                        if a.ISYMBL not in imp]
        watoms.sort(key=lambda a: a.GetIdx())

        atoms = watoms + eatoms

    # Put back the central atom
    atoms.insert(2, c)

    # Check connectivity
    if (not atoms[0].IsConnected(atoms[1])
        and atoms[2].IsConnected(atoms[0])
        and atoms[2].IsConnected(atoms[1])
        and atoms[2].IsConnected(atoms[3])):
        return atoms
    else:
        return None


# TODO: Add an option to set improper from PREP file even if it doesn't match
def findImpropers(mol, imps):
    """Find impropers in molecule.

    """
    smarts = ("[CD3](~[OD1])(~[OD1])~[*]", #  X-O2- C-O2
              "[ND3](~[HD1])(~[*])~[*]",   #  X- X- N- H
              "[ND3](~[CD4])(~[*])~[*]",   #  X- X- N-CG
              "[CD3](~[OD1])(~[*])~[*]",   #  X- X- C- O
              "[CD3](~[OD1])(~[OD2])~[*]", #  X-O2- C-OH
              "[ND3](~[HD1])(~[S])~[CD4]", # CG- S- N- H
              "[CD3](~[*])(~[*])~[*]",     # sp2 C
             )
    impropers = []
    sm = ob.OBSmartsPattern()

    for i in smarts:
        if not sm.Init(i):
            raise IOError("Invalid SMARTS pattern")
        sm.Match(mol)
        matches = list(sm.GetUMapList()) # List of atom numbers
        # Check each of them for matching atomtypes in imps
        for j in matches:
            found = False
            types = {mol.GetAtom(a).ISYMBL for a in j}
            for imp in imps.keys():
                # If equivalent or
                # if all atomtypes from imp are present
                impSet = set(imp)
                impSet.discard('X')
                if impSet.issubset(types):
                    # Check bonding and reorder atoms
                    atoms = getOrderedImp([mol.GetAtom(a) for a in j], imp)
                    if atoms is not None:
                        # Add to impropers
                        impropers.append([a.GetIdx() for a in atoms] +
                                         [imps[imp][0]])
                        # Move to next match
                        found = True
                        break
            if not found:
                # It was not in forcefield, threat it as normal Gromacs imp
                imp = list(j)[1:]
                imp.sort()
                imp.insert(0, j[0])
                impropers.append(imp + [2])

    impropers = uniquify(impropers)
    impropers.sort()
    if impropers != []:
        print("Impropers found: {}".format(len(impropers)))

    return impropers


def getV2(sigma1, sigma2):
    """Combine sigmas using combination rule 2.

    """
    return (sigma1 + sigma2) / 2


def getW2(epsilon1, epsilon2, scnb=1):
    """Combine epsilons using combination rule 2.

    """
    return sqrt(epsilon1 * epsilon2) / scnb


def getBonds(mol):
    """Get list of bonds in the molecule.

    """
    bonds = []
    for bond in ob.OBMolBondIter(mol):
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        s = [a, b]
        s.sort()
        bonds.append(s)
    bonds.sort()

    return bonds


def getAngles(mol):
    """Get list of angles in the molecule.

    """
    angles = []
    for angle in ob.OBMolAngleIter(mol):
        angles.append([angle[1]+1, angle[0]+1, angle[2]+1])
    for i in angles:
        if i[2] < i[0]:
            i.reverse()
    angles.sort(key=lambda a: [a[1], a[0], a[2]])

    return angles


def getTorsions(mol):
    """Get list of torsions in the molecule.

    """
    torsions = []
    for torsion in ob.OBMolTorsionIter(mol):
        torsions.append([torsion[0]+1, torsion[1]+1,
                         torsion[2]+1, torsion[3]+1])
    for i in torsions:
        if i[2] < i[1]:
            i.reverse()
    torsions.sort(key=lambda a: [a[1], a[2], a[0], a[3]])

    return torsions


def getPairs(torsions, angles, bonds):
    """Get list of 1-4 pairs in the molecule.

    """
    pairs14 = [sorted([i[0], i[3]]) for i in torsions]
    # Filter out duplicate pairs and those in 4 and 5 member rings
    excludel = [[j[0],j[2]] for j in angles]
    excludel.extend(bonds)
    pairs14 = uniquify(list_difference(pairs14, excludel))
    pairs14.sort()

    return pairs14


def calc14NBParams(mol, pairs14, sc14):
    """Calculate explicit 14-LJ parameters for function type 2.

    """
    aNames = {a.GetIdx(): a.ISYMBL for a in ob.OBMolAtomIter(mol)}
    for pair in pairs14:
        for d in sc14[0].keys():
            if ({aNames[pair[0]], aNames[pair[1]]} == {d[0], d[3]}):
                # Found
                try:
                    scee, scnb = sc14[0][d]
                    s1, e1 = sc14[1][d[0]]
                    s2, e2 = sc14[1][d[3]]
                except KeyError:
                    # Some of the atoms didn't have LJ parameters in the
                    # Glycam ITP file, ignore them
                    break
                q1 = mol.GetAtom(pair[0]).GetPartialCharge()
                q2 = mol.GetAtom(pair[1]).GetPartialCharge()
                V = getV2(s1, s2)
                W = getW2(e1, e2, scnb)

                # Calculate explicit parameters
                pair.extend([1/scee, q1, q2, V, W])

                # Move to the next pair
                break


def getExplicitImps(mol):
    """Get list of explicit impropers from PREP file(s).

    """
    impropers = []
    hasImpropers = False
    for i in ob.OBResidueIter(mol):
        if i.HasData('IMPROPER'):
            hasImpropers = True
            imps2 = i.GetData('IMPROPER').GetValue().split()
            for j in imps2:
                imp = []
                atoms = j.split('-')
                for a in atoms:
                    for k in ob.OBResidueAtomIter(i):
                        if k.IGRAPH == a:
                            imp.append(k.GetIdx())
                if len(imp) != 4:
                    s = "Error: Could not find improper {0} in residue {1}"
                    print(s.format(j, i.GetName()))
                # To match inconsistency in Leap
                if atoms[1][0] == 'N':
                    imp[0], imp[1] = imp[1], imp[0]
                impropers.append(imp + [4])
    impropers.sort()

    return impropers, hasImpropers

def suffixSugars(mol, n1, n2):
    """Add suffix '_S' to atomtypes for all atoms in given range.

    """
    s1 = set(range(n1, n2))
    for atom in ob.OBMolAtomIter(mol):
        if atom.GetIdx() in s1:
            atom.ISYMBL += "_S"


def writeGromacsTop (mol, itpFN, impTypes={}, sc14={}):
    """Write Gromacs topology.

    This function takes an OBMol object and writes the corresponding gromacs
    topology file (itp).
    """
    nexcl = 3 # minim distance between atoms to have long range interactions
    # Bonded parameter function types
    bfunc = 1 # bond
    pfunc = 1 # 1-4 pair
    afunc = 1 # angle
    dfunc = 9 # dihedral

    bonds    = getBonds(mol)
    angles   = getAngles(mol)
    torsions = getTorsions(mol)

    pairs14  = getPairs(torsions, angles, bonds)

    calc14NBParams(mol, pairs14, sc14)

    # Find impropers in molecule
    imps = findImpropers(mol, impTypes)

    # Get explicit impropers from PREP file if any
    impropers, hasImpropers = getExplicitImps(mol)

    # Check if all explicit impropers were found
    if (hasImpropers
        and not set(map(tuple, impropers)).issubset(set(map(tuple, imps)))):
        print("Error: Did not find all explicit impropers!")
        for i in imps:
            print("found:", i)
        for i in impropers:
            print("explicit:", i)
        sys.exit(2)

    # Use found impropers
    if imps != []:
        impropers = imps
        hasImpropers = True


    if (options.opls):
        if options.amber:
            options.opls == False
        else:
            print(" Suffixing polymer atomtypes ...")
            suffixSugars(mol, 1, mol.NumAtoms()+1)
            print(" Done.")

    # Construct the itp file
    t = []
    t.append("[ moleculetype ]\n")
    t.append("; Name    nrexcl\n")
    props = {"MolNM": mol.GetTitle(), "nrexcl": 3}
    t.append("{0[MolNM]:s} {0[nrexcl]:d}\n".format(props))

    t.append("\n[ atoms ]\n")
    t.append(";   nr       type      resnr   residue      atom  cgnr    charge"
             "     mass\n")
    s = "{0[aID]:6d}{0[aT]:>11s}{0[resID]:7d}{0[resN]:>7s}" \
        "{0[aN]:>7s}{0[aCG]:7d}{0[aQ]:11.4f}{0[aM]:11.3f}\n"
    for i in ob.OBResidueIter(mol):
        t.append("; Residue {0:d} {1:s} q={2:.3f}\n".format(i.GetIdx()+1,
                                                            i.GetName(),
                                                            getResCharge(i)))
        t.append("; {0:s}\n".format(i.GetData("NOTE").GetValue()))
        props = {"resID": i.GetIdx()+1, "resN": i.GetName()}
        for j in ob.OBResidueAtomIter(i):
            props["aID"] = j.GetIdx()
            props["aT"]  = j.ISYMBL
            props["aN"]  = j.IGRAPH
            props["aCG"] = j.GetIdx()
            props["aQ"]  = j.GetPartialCharge()
            props["aM"]  = j.GetAtomicMass()
            t.append(s.format(props))

    # Bonds
    t.append("\n[ bonds ]\n")
    t.append(";   ai     aj   funct   c0   c1\n")
    for i in bonds:
        t.append("{0[0]:6d} {0[1]:6d} {1:6d}\n".format(i, bfunc))

    # 1-4 pairs
    t.append("\n[ pairs ]\n")
    t.append(";   ai     aj   funct   c0   c1\n")
    s = "{0[0]:6d} {0[1]:6d} {1:6d}   {2:6.4f} {3:9.4f} {4:9.4f}" \
        "   {5:1.5e}   {6:1.5e}\n"
    for i in pairs14:
        if len(i) > 2:
            t.append(s.format(i, 2, i[2], i[3], i[4], i[5], i[6]))
        else:
            t.append("{0[0]:6d} {0[1]:6d} {1:6d}\n".format(i, pfunc))

    # Angles
    t.append("\n[ angles ]\n")
    t.append(";    ai    aj     ak   funct   c0   c1\n")
    for i in angles:
        t.append("{0[0]:6d} {0[1]:6d} {0[2]:6d} {1:6d}\n".format(i, afunc))

    # Proper torsions
    t.append("\n[ dihedrals ]\n")
    t.append(";   ai     aj     ak     al  funct\n")
    s = "{0[0]:6d} {0[1]:6d} {0[2]:6d} {0[3]:6d} {1:6d}\n"
    for i in torsions:
        t.append(s.format(i, dfunc))

    # Improper torsions
    if hasImpropers:
        t.append("\n[ dihedrals ]\n")
        t.append(";   ai     aj     ak     al  funct\n")
        for i in impropers:
            t.append(s.format(i, i[4]))

    t.append("\n")

    # Write the file
    with open(itpFN, 'w') as f:
        f.writelines(t)


def distortMol(mol, d, n):
    conf = []
    for i in range(n):
        coord = []
        for a in ob.OBMolAtomIter(mol):
            coord.append([a.x()+random.uniform(-d, d),
                          a.y()+random.uniform(-d, d),
                          a.z()+random.uniform(-d, d)])
        conf.append(coord)
    return conf


def generateMonomers(residues,
                     fformat='pdb',
                     impTypes={},
                     sc14={},
                     distortion=0.0):
    """Generate PDBs for all fragments whose name starts with 0.

    """
    skipRes = ["0DD", "0GL"] # Missing some parameters in Amber
    for k in residues.keys():
        if k in skipRes:
            continue
        if k[0] == '0':
            at = residues[k].GetFirstAtom().IGRAPH
            mol = buildSeq(residues, ["{0}=ROH-(O1,{1})-{0}".format(k, at)])
            mol.SetTitle("{0}".format(k))
            writeMol(mol, fformat, k + '_000.' + fformat)
            writeGromacsTop(mol, "{0}_{1:03d}.itp".format(k, 0),
                            impTypes=impTypes, sc14=sc14)
            if distortion != 0.0:
                confs = distortMol(ob.OBMol(mol), distortion, 10)
                for n, i in enumerate(confs):
                    for a in ob.OBMolAtomIter(mol):
                        a.SetVector(i[a.GetIdx()-1][0],
                                    i[a.GetIdx()-1][1],
                                    i[a.GetIdx()-1][2])
                    writeMol(mol, fformat,
                             "{0}_{1:03d}.{2}".format(k, n+1, fformat))
                    writeGromacsTop(mol, "{0}_{1:03d}.itp".format(k, n+1),
                                    impTypes=impTypes, sc14=sc14)


def generateDimers(residues,
                   fformat='pdb',
                   impTypes={},
                   sc14={},
                   distortion=0.0):
    """Generate PDBs and ITPs for all homodimers.

    """
    skipRes = ["2SA"] # OB can't generate 3D structure for 06g
    s = "{2}=ROH-(O1,{1})-{2}-(O{2[0]},{3})-{0}"
    for k in residues.keys():
        if k[0] == '0':
            a1 = residues[k].GetFirstAtom().IGRAPH
            for i in (1,2,3,4,6):
                d = str(i) + k[1:]
                if d in skipRes:
                    continue
                try:
                    a2 = residues[d].GetFirstAtom().IGRAPH
                except KeyError:
                    continue
                mol = buildSeq(residues, [s.format(k, a2, d, a1)])
                mol.SetTitle("{}".format(d))
                writeMol(mol, fformat, "{0}_{1:03d}.{2}".format(d, 0, fformat))
                writeGromacsTop(mol, "{0}_{1:03d}.itp".format(d, 0),
                                impTypes=impTypes, sc14=sc14)
                if distortion != 0.0:
                    confs = distortMol(ob.OBMol(mol), distortion, 10)
                    for n, i in enumerate(confs):
                        for a in ob.OBMolAtomIter(mol):
                            a.SetVector(i[a.GetIdx()-1][0],
                                        i[a.GetIdx()-1][1],
                                        i[a.GetIdx()-1][2])
                        writeMol(mol, fformat,
                                 "{0}_{1:03d}.{2}".format(d, n+1, fformat))
                        writeGromacsTop(mol, "{0}_{1:03d}.itp".format(d, n+1),
                                        impTypes=impTypes, sc14=sc14)


def generateTrimers(residues,
                    fformat='pdb',
                    impTypes={},
                    sc14={},
                    distortion=0.0):
    """Generate PDBs and ITPs for all homotrimers.

    """
    skipRes = ["2SA"] # OB can't generate 3D structure for 06g
    s = "{2}{4}=ROH-(O1,{1})-{2}-(O{2[0]},{3})-{4}-(O{4[0]},{5})-{0}"
    for k in residues.keys():
        if k[0] == '0':
            a1 = residues[k].GetFirstAtom().IGRAPH
            for i in (1,2,3,4,6):
                d = str(i) + k[1:]
                if d in skipRes:
                    continue
                try:
                    a2 = residues[d].GetFirstAtom().IGRAPH
                except KeyError:
                    continue
                for j in (1,2,3,4,6):
                    t = str(j) + k[1:]
                    if t in skipRes:
                        continue
                    try:
                        a3 = residues[t].GetFirstAtom().IGRAPH
                    except KeyError:
                        continue
                    mol = buildSeq(residues,
                                   [s.format(k, a2, d, a3, t, a1)])
                    mol.SetTitle("{}{}".format(d, t))
                    writeMol(mol, fformat,
                             "{0}_{1:03d}.{2}".format(d+t, 0, fformat))
                    writeGromacsTop(mol, "{0}_{1:03d}.itp".format(d+t, 0),
                                    impTypes=impTypes, sc14=sc14)
                    if distortion != 0.0:
                        confs = distortMol(ob.OBMol(mol), distortion, 10)
                        for n, i in enumerate(confs):
                            for a in ob.OBMolAtomIter(mol):
                                a.SetVector(i[a.GetIdx()-1][0],
                                            i[a.GetIdx()-1][1],
                                            i[a.GetIdx()-1][2])
                            writeMol(mol, fformat,
                                     "{0}_{1:03d}.{2}".format(d+t, n+1,
                                                              fformat))
                            writeGromacsTop(mol,
                                            "{0}_{1:03d}.itp".format(d+t, n+1),
                                            impTypes=impTypes, sc14=sc14)


def writeMol(mol, fformat, fname):
    """Write coordinate file for a given residue.

    """
    converter = ob.OBConversion()
    converter.SetOutFormat(fformat)
    converter.WriteFile(mol, fname)


def charmmGui (inputPDB):
    atomPDB = []
    hetatmPDB = []
    conect = []
    header = "TITLE CHARMM-GUI COMPATIBLE\n"
    boxsize = []
    box = re.compile(r"CRYST")

    with open(inputPDB) as file:    
        for line in file:
            if line.startswith("ATOM"):
                atomPDB.append(line)
            if line.startswith("HETATM"):
                hetatmPDB.append(line)
            if line.startswith("CONECT"):
                conect.append(line)
            if re.match(box,line):
                    boxsize.append(line)
    
    
    # Replace ATOM with HETATM and viceversa:
    hetatm      = [w.replace("ATOM  " , "HETATM") for w in atomPDB]
    #new_ATOM            = [w.replace("HETATM" , "ATOM  ") for w in HETATM]
        
    output = str(inputPDB).strip(".pdb") + "_charmm-gui.pdb"
    
    with open(output,"w") as newfile:    
        newfile.write(header)
        newfile.write("".join(boxsize))
        newfile.write("".join(hetatm))
        newfile.write("".join(conect))
        newfile.write("END")



if __name__ == "__main__":
    options, args = optP()

    if (options.ffFN is None and
        (options.itpFN is not None
         or options.gen_mono
         or options.gen_di
         or options.gen_tri)):
        print("Error: Glycam forcefield ITP file has to be specified!")
        sys.exit(1)

    if options.ffFN is not None and not os.path.isfile(options.ffFN):
        print("Error: File {} not found!".format(options.ffFN))
        sys.exit(1)

    if (options.molFN is None
        and not options.gen_mono
        and not options.gen_di
        and not options.gen_tri):
        print("Error: Specify filename for structure!")
        sys.exit(1)

    if options.ffFN is not None:
        sc14 = get14NBParams(options.ffFN)
        imps = getImproperTypes(options.ffFN)

    residues = {}
    for i in options.prepFNs:
        print(i)
        with open(i, 'r') as f:
            ires = readPrepFile(f, options.no_warn)
        for k, v in ires.items():
            residues[k] = v

    print("Got", len(residues), "molecules")

    if (options.gen_mono):
        generateMonomers(residues, impTypes=imps, sc14=sc14, distortion=0.2)
        print("Done.")
        sys.exit(0)

    if options.gen_di:
        generateDimers(residues, impTypes=imps, sc14=sc14, distortion=0.2)
        print("Done.")
        sys.exit(0)

    if options.gen_tri:
        generateTrimers(residues, impTypes=imps, sc14=sc14, distortion=0.2)
        print("Done.")
        sys.exit(0)

    with open(options.seqFN, 'r') as f:
        sequence = []
        lines = f.readlines()
        for i in lines:
            sequence.append(i.rstrip('\n'))
    mol = buildSeq(residues, sequence)
    filename = options.molFN
    fformat = os.path.splitext(filename)[1][1:]
    writeMol(mol, fformat, filename)

    if (options.itpFN):
        writeGromacsTop(mol, options.itpFN, impTypes=imps, sc14=sc14)

    if (options.pathFN):
        print ("Writing Charmm-Gui compatible PDB file")
        if  fformat != "pdb":
            print ("\nWARNING: your output File has to be in PDB format!\n"
                   "The Charmm-Gui compatible PDB file has not been created") 
        else :
            charmmGui(options.molFN)
            print ("\nNOTE : The Charmm-Gui won't read structures with steric clashes.\n" \
               "Please, consider to minimize your output structure first and then run\n" \
               "doglycans2charmmgui.py separately.")
            print ("Done.")