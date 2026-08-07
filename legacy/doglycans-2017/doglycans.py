#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Version 1
Authors: See the AUTHORS file for a list of contributors to doGlycans.py
License GPLv3

Attach sugar chains to protein.

Right now there are still quite a lot hardcoded assumptions about the system:
    - Protein chain in sequence files is referenced as "P"
    - All new atoms come after the protein atoms in topology
    - Impropers for sugar and cross-link part have only generic values
    - The bonding in protein is inferred from atom distances
    - Conformation of sugar chains is read from sequence files

Input:
    - Protein structure (PDB)
    - Protein topology (for reading impropers and verifying bonding)
    - Sequence file for sugars with specified attachment points to protein
    - Glycam PREP file(s) with sugar residues
    - Glycam forcefield ITP file

    + Mark sugar cross-links with protein
    + Attach impropers from topology to residues
    + Mark sugar only pairs for 1-4 scaling correction

Output:
    - Glycosylated protein structure
    - Gromacs topology file with marked protein-sugar cross-link parts

Manual modifications still required:
    - Changing atom types in protein for attachment point N and H
    - Changing charges for N and H (+ adjust to get integer charge for protein)
    - Adding explicit bonded parameters for cross-link
"""

import re
import sys
import os.path
import pybel
import numpy as np
import openbabel as ob
import prepreader as pr

from collections import OrderedDict
from optparse import OptionParser


def optP():
    """Parse command line options.

    """
    def getFNs(option, opt_str, value, parser):
            """Extract numbers from command line.

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

    usage="[python3] %prog [-W] -f $DGPATH/dat/glycam06h.itp" \
          " -p $DGPATH/dat/GLYCAM_06h.prep -s sequence_file.seq" \
          " -c ProtH.pdb -i ProtH.itp -o ProtG.pdb -t ProtG.itp -B" \
          "  --opls --amber --charmm $DGPATH/dat/GLYCAM_06h_resname.txt"

    description="Attach sugar chains to protein."
    version="\n%prog Version 0.04\n\nRequires Python 3, Open Babel Python" \
            " bindings, numpy and prepReader.py."

    optParser = OptionParser(usage=usage,
                             version=version,
                             description=description)

    optParser.set_defaults(ffFN=None, prepFN=['Glycam_06.prep'], seqFN=None,
                           protFN=None, postFN=None, itpFNin=None, itpFNout=None,
                           noBonds=False, prefix=False, opls=True, charmm=None,
                           amber=False, no_warn=False)

    optParser.add_option('-f', type='str',
                         dest='ffFN',
                         help="Glycam parameter ITP file (in)"
                         " [default: %default]")
    optParser.add_option('-p', action='callback', callback=getFNs,
                         dest='prepFNs',
                         help="Prep file(s) [default: %default]")
    optParser.add_option('-s', type='str',
                         dest='seqFN',
                         help="Sequence file [default: %default]")
    optParser.add_option('-c', type='str',
                         dest='protFN',
                         help="Protein structure file [default: %default]")
    optParser.add_option('-o', type='str',
                         dest='postFN',
                         help="Output structure file [default: %default]")
    optParser.add_option('-W', action='store_true', default=False,
                         dest='no_warn',
                         help="Suppress warnings [default: %default]")
    optParser.add_option('-i', type='str',
                         dest='itpFNin',
                         help="Gromacs topology file (in) [default: %default]")
    optParser.add_option('-t', type='str',
                         dest='itpFNout',
                         help="Gromacs topology file (out) [default: %default]")
    optParser.add_option('-g', action='store_true',
                         dest='prefix',
                         help="Add prefix 'G_' for atom types"
                         " [default: %default]")
    optParser.add_option('-B', action='store_true',
                         dest='noBonds',
                         help="Don't perceive bonds from atom distances"
                         " [default: %default]")
    optParser.add_option('--opls', action='store_true',
                         dest='opls',
                         help="Add suffix '_S' for opls atom types"
                         " [default: %default]")
    optParser.add_option('--charmm', type='str',
                         dest='pathFN',
                         help="Print PDB for charmm-gui submission.\n" \
                         "The path to a file with sugar residue name list is needed"
                         " [default: %default]")
    optParser.add_option('--amber', action='store_true',
                         dest='amber',
                         help="Remove suffix '_S' for amber atom types"
                         " [default: %default]")

    options, args = optParser.parse_args()

    if options.ffFN is None:
        s = "Error: Glycam ITP file required."
        print(s)
        sys.exit(2)

    if options.protFN is None:
        s = "Error: Protein structure file required."
        print(s)
        sys.exit(2)

    if options.postFN is None:
        s = "Error: Output structure file required."
        print(s)
        sys.exit(2)

    if options.seqFN is None:
        s = "Error: Sugar sequence file required."
        print(s)
        sys.exit(2)

    if options.itpFNin is None:
        s = "Error: Input ITP file required."
        print(s)
        sys.exit(2)

    if options.itpFNout is None:
        s = "Error: Output ITP file name required."
        print(s)
        sys.exit(2)

    return options, args


def turnUp(mol, a, b, c, d):
    """Rotate dihedral to point bond a-b away from the geometric center.

    """
    # Get the geometric center
    x = 0.0
    y = 0.0
    z = 0.0
    for atom in ob.OBMolAtomIter(mol):
        x += atom.x()
        y += atom.y()
        z += atom.z()
    x /= mol.NumAtoms()
    y /= mol.NumAtoms()
    z /= mol.NumAtoms()
    CM = np.array([x, y, z])

    # Calculate the angle
    ac = np.array([a.x(), a.y(), a.z()])
    bc = np.array([b.x(), b.y(), b.z()])
    cc = np.array([c.x(), c.y(), c.z()])

    r1 = cc - CM
    r2 = ac - bc
    r3 = bc - cc

    nr1 = np.linalg.norm(r1)
    nr2 = np.linalg.norm(r2)
    nr3 = np.linalg.norm(r3)

    co    = np.dot(r1, r2)/(nr1*nr2)
    angle = np.arccos(np.clip(co, -1, 1))
    print("cos θ =", co)
    print("θ     = {: .1f}".format(np.degrees(angle)))

    si    = np.linalg.norm(np.cross(r1, r2))/(nr1*nr2)
    angle = np.arcsin(np.clip(si, -1, 1))
    sign  = np.dot(np.cross(r1, r2)/(nr1*nr2), r3/nr3)
    print("sin θ =", si)
    print("θ     = {: .1f}".format(np.degrees(angle)))
    print("Sign:", sign)

    if abs(si) < 1e-14 and abs(co) < 1e-14:
        angle = 0.0
    elif sign < 0:
        angle = -np.arctan2(si, co)
    else:
        angle = np.arctan2(si, co)
    print("θ     = {: .1f}".format(np.degrees(angle)))

    # Get the torsion value
    dih = np.radians(mol.GetTorsion(d, c, b, a))
    print("Torsion: {: .1f}".format(mol.GetTorsion(d, c, b, a)))
    print("Rotate:  {: .1f}".format(np.degrees(angle)))

    # Set the torsion
    dih -= angle
    print("Rotate dihedral in residue by {:.1f} deg".format(np.degrees(angle)))
    mol.SetTorsion(d, c, b, a, dih)
    print("Torsion: {: .1f}".format(mol.GetTorsion(d, c, b, a)))


def buildSeq(residues, parent, sequence):
    """Build a molecule with given sequence on parent molecule.

    """
    builder = ob.OBBuilder()
    # Use ring conformations from existing 3D coordinates
    builder.SetKeepRings()
    parent.SetAutomaticPartialCharge(False)

    print("Seq:")
    for i in sequence:
        print(i)

    # TODO: Change atom type for N and H at linking?
    for nres in pr.getSeqRes(sequence):
        if nres.name is not None:
            # Add next residue
            parent += residues[nres.name]

            # Have to manually copy residue data, because += operator doesn't
            res = parent.GetResidue(parent.NumResidues()-1)
            # Increment residue number from the previous residue
            res.SetNum(parent.GetResidue(parent.NumResidues()-2).GetNum()+1)
            res.CloneData(residues[nres.name].GetResidue(0).GetData('NOTE'))
            res.CloneData(residues[nres.name].GetData('IMPROPER'))
            pr.addPairData(res, 'TRUNK_ID', nres.trunk_id)
            pr.addPairData(res, 'BRANCH_ID', nres.branch_id)

            # Get parent residue from the right chain
            p = pr.getResByBranchID(parent, nres.trunk_id)

            # Remove the atoms replaced by attachment if any
            if nres.replace is not None:
                if nres.replace == 'HD22':
                    a = pr.getResAtomByName(p, 'ND2')
                    b = pr.getResAtomByName(p, 'CG')
                    c = pr.getResAtomByName(p, 'CB')
                    d = pr.getResAtomByName(p, 'CA')
                    turnUp(parent, a, b, c, d)
                    a = pr.getResAtomByName(p, 'HD22')
                    b = pr.getResAtomByName(p, 'ND2')
                    c = pr.getResAtomByName(p, 'CG')
                    d = pr.getResAtomByName(p, 'CB')
                    turnUp(parent, a, b, c, d)
                atom = pr.getResAtomByName(p, nres.replace)
                print("Replaced atom {} {}".format(atom.GetIdx(),
                                                   nres.replace))
                pr.addPairData(p,
                               'REPLACED_ATOM',
                               "{}:{}".format(p.GetNum(), nres.replace))
                p.RemoveAtom(atom)
                parent.DeleteAtom(atom)

            # Get indices for atoms in connection
            a = pr.getResAtomByName(p,   nres.bond[0]).GetIdx()
            b = pr.getResAtomByName(res, nres.bond[1]).GetIdx()

            # Join the residues
            if builder.Connect(mol, a, b):
                # Set torsion(s)
                pr.setTorsions(mol, p, nres)
            else:
                print("Error: Something went wrong while joining fragments!")
                sys.exit(3)
        else:
            # Don't add new residue, only make given connection
            r = pr.getResByBranchID(mol, nres.branch_id)

            # Get parent residue
            p = pr.getResByBranchID(mol, nres.trunk_id)

            # Get atoms
            a = pr.getResAtomByName(p, nres.bond[0])
            b = pr.getResAtomByName(r, nres.bond[1])

            # Make the bond
            pr.addBond(mol, a, b)

            # Set torsion(s)
            pr.setTorsions(mol, p, nres)

    return parent


def getAtomTypes(f):
    """Get atom types from topology.

    """
    reXa = re.compile(r'^[^;]*\[\s+atoms\s+\]')
    reXc = re.compile(r'^\s*[;]+')
    reXe = re.compile(r'^\s*$')
    reXs = re.compile(r'^[^;]*\[\s+[a-z]*\s+\]')
    atypes     = {}
    resnumbers = {}
    resnames   = {}
    anames     = {}
    chgroups   = {}
    charges    = {}
    line = f.readline()
    while reXa.match(line) is None:
        line = f.readline()

    line = f.readline()
    while reXs.match(line) is None:
        if (reXc.match(line) is not None
            or reXe.match(line) is not None):
            line = f.readline()
            continue
        anr, atype, resnr, resn, aname, chgr, ch, mass = line.split()[:8]
        atypes[int(anr)]     = atype
        resnumbers[int(anr)] = int(resnr)
        resnames[int(resnr)] = resn
        anames[int(anr)]     = aname
        chgroups[int(anr)]   = chgr
        charges[int(anr)]    = ch
        line = f.readline()

    return atypes, resnumbers, resnames, anames, chgroups, charges


def getFunctionTypes(f):
    """Get function types from topology.

    """
    reXb = re.compile(r'^[^;]*\[\s+bonds\s+\]')
    reXp = re.compile(r'^[^;]*\[\s+pairs\s+\]')
    reXa = re.compile(r'^[^;]*\[\s+angles\s+\]')
    reXd = re.compile(r'^[^;]*\[\s+dihedrals\s+\]')
    reXc = re.compile(r'^\s*[;]+')
    reXe = re.compile(r'^\s*$')
    reXs = re.compile(r'^[^;]*\[\s+[a-z]*\s+\]')

    bfn = None # bonds
    pfn = None # pairs
    afn = None # angles
    dfn = None # dihedrals
    ifn = None # impropers

    f.seek(0)
    line = f.readline()
    while (reXb.match(line) is None
           and line != ''):
        line = f.readline()

    line = f.readline()
    while (reXs.match(line) is None
           and line != ''):
        if (reXc.match(line) is not None
            or reXe.match(line) is not None):
            line = f.readline()
            continue
        bfn = int(line.split()[2])
        break

    f.seek(0)
    line = f.readline()
    while (reXp.match(line) is None
           and line != ''):
        line = f.readline()

    line = f.readline()
    while (reXs.match(line) is None
           and line != ''):
        if (reXc.match(line) is not None
            or reXe.match(line) is not None):
            line = f.readline()
            continue
        pfn = int(line.split()[2])
        break

    f.seek(0)
    line = f.readline()
    while (reXa.match(line) is None
           and line != ''):
        line = f.readline()

    line = f.readline()
    while (reXs.match(line) is None
           and line != ''):
        if (reXc.match(line) is not None
            or reXe.match(line) is not None):
            line = f.readline()
            continue
        afn = int(line.split()[3])
        break

    f.seek(0)
    line = f.readline()
    while (reXd.match(line) is None
           and line != ''):
        line = f.readline()

    line = f.readline()
    while (reXs.match(line) is None
           and line != ''):
        if (reXc.match(line) is not None
            or reXe.match(line) is not None):
            line = f.readline()
            continue
        dfn = int(line.split()[4])
        break

    line = f.readline()
    while (reXd.match(line) is None
           and line != ''):
        line = f.readline()

    line = f.readline()
    while (reXs.match(line) is None
           and line != ''):
        if (reXc.match(line) is not None
            or reXe.match(line) is not None):
            line = f.readline()
            continue
        ifn = int(line.split()[4])
        break

    return bfn, pfn, afn, dfn, ifn


def getITPBonds(f):
    """Get bond list from topology.

    """
    reXb = re.compile(r'^[^;]*\[\s+bonds\s+\]')
    reXc = re.compile(r'^\s*[;]+')
    reXe = re.compile(r'^\s*$')
    reXs = re.compile(r'^[^;]*\[\s+[a-z]*\s+\]')

    bonds = []

    f.seek(0)
    line = f.readline()
    while (reXb.match(line) is None
           and line != ''):
        line = f.readline()

    line = f.readline()
    while (reXs.match(line) is None
           and line != ''):
        if (reXc.match(line) is not None
            or reXe.match(line) is not None):
            line = f.readline()
            continue
        a, b = sorted(line.split()[:2])
        bonds.append([int(a), int(b)])
        line = f.readline()
    bonds.sort()

    return bonds


def getImpropers(f):
    """Get impropers from topology.

    """
    reXa = re.compile(r'^[^;]*\[\s+dihedrals\s+\]')
    reXc = re.compile(r'^\s*[;]+')
    reXe = re.compile(r'^\s*$')
    line = f.readline()

    impropers = []
    improper_atoms = set()
    while line != '':
        while reXa.match(line) is None:
            line = f.readline()
            if line == '':
                break
        if line == '':
            break
        line = f.readline()

        while reXe.match(line) is None:
            if reXc.match(line) is not None:
                line = f.readline()
                continue
            a, b, c, d, func = line.split()[:5]
            if func == '1':
                impropers.append([int(a),
                                  int(b),
                                  int(c),
                                  int(d),
                                  int(func),
                                  line.split()[5:]])
                for i in impropers[-1][:4]:
                    improper_atoms.add(i)
            if func == '4':
                impropers.append([int(a),
                                  int(b),
                                  int(c),
                                  int(d),
                                  int(func)])
                for i in impropers[-1][:4]:
                    improper_atoms.add(i)
            line = f.readline()

    return impropers, improper_atoms


def getDefTorsions(f):
    """Get defined dihedrals from topology.

    """
    reXd = re.compile(r'^[^;]*\[\s+dihedrals\s+\]')
    reXp = re.compile(r'^torsion_[A-Z]{3}_[A-Z0-9_]+_mult[0-9]')
    reXc = re.compile(r'^\s*[;]+')
    reXe = re.compile(r'^\s*$')
    line = f.readline()

    defTorsions = {}
    while line != '':
        while (reXd.match(line) is None
               and line != ''):
            line = f.readline()
        if line == '':
            break
        line = f.readline()

        while reXe.match(line) is None:
            if (reXc.match(line) is not None
                or reXe.match(line) is not None):
                line = f.readline()
                continue
            if len(line.split()) > 5:
                a, b, c, d, func, rem = line.split()[:6]
                a, b, c, d = int(a), int(b), int(c), int(d)
                if reXp.match(rem) is not None:
                    if (a,b,c,d) in defTorsions:
                        defTorsions[(a,b,c,d)].append(rem)
                    else:
                        defTorsions[(a,b,c,d)] = [rem]
            line = f.readline()

    return defTorsions


def setDefTorsions(torsions, defTorsions):
    """Add defined torsions to the dihedrals list.

    """
    s = "{}:{}"
    deleteDihs = []
    tmpTorsions = []
    for torsion in torsions:
        atom = mol.GetAtom(torsion[0])
        a = s.format(atom.GetResidue().GetNum(), atom.IGRAPH)
        atom = mol.GetAtom(torsion[1])
        b = s.format(atom.GetResidue().GetNum(), atom.IGRAPH)
        atom = mol.GetAtom(torsion[2])
        c = s.format(atom.GetResidue().GetNum(), atom.IGRAPH)
        atom = mol.GetAtom(torsion[3])
        d = s.format(atom.GetResidue().GetNum(), atom.IGRAPH)
        if (a,b,c,d) in defTorsions:
            t = [torsion[:4] for i in range(len(defTorsions[(a,b,c,d)]))]
            deleteDihs.append(torsion)
            for n, i in enumerate(defTorsions[(a,b,c,d)]):
                t[n].append(i)
            tmpTorsions.extend(t)
    torsions.extend(tmpTorsions)
    for i in deleteDihs:
        torsions.remove(i)
    torsions.sort(key=lambda d: [d[1], d[2], d[0], d[3], d[4]])


def mark_cross(mol, s1, s2, l):
    """Mark bonding lists for crossection between protein and sugar.

    """
    s3 = set()
    remf = " {0:>9s}"
    for i in l:
        si = set(i)
        if (not si.isdisjoint(s1)
            and not si.isdisjoint(s2)):
            s3.add(str(i))
    for k, j in enumerate(l):
        if str(j) in s3:
            rems = " ; mark"
            for i in l[k]:
                a =  mol.GetAtom(i)
                ra = a.GetResidue()
                rems += remf.format(str(ra.GetNum()) + ':' + ra.GetAtomID(a))
            l[k].append(rems)
        else:
            l[k].append('')


def getResidueImps(mol, impropers, hasImpropers):
    """Add impropers stored in residues (protein and marked).

    """
    marked_impropers = []
    for res in ob.OBResidueIter(mol):
        if res.HasData('IMPROPERS'):
            hasImpropers = True
            imps = res.GetData('IMPROPERS').GetValue().split('-')
            if imps == ['']:
                s = "Error: Empty impropers in residue {}"
                print(s.format(res.GetNum()))
                sys.exit(3)
            for i in imps:
                imp = []
                natoms = i.split()[:4]
                # Check if the atom has been removed
                if (res.HasData('REPLACED_ATOM') and
                    res.GetData('REPLACED_ATOM').GetValue() in natoms):
                    # Get the name of replaced atom
                    ra = res.GetData('REPLACED_ATOM').GetValue()
                    rar, ran = ra.split(':')
                    for k, n in enumerate(natoms):
                        if n == ra:
                            # Get the index of replaced atom in the improper
                            nr = k
                        else:
                            # Get the neighbor of replaced atom
                            nbr = pr.getResAtomByName(res,
                                                      n.split(':')[1])
                    rn = pr.getAtomNbrByName(nbr, "C1").GetResidue().GetNum()
                    natoms[nr] = str(rn) + ":C1"

                    # Atom indices for the improper
                    nnatoms = []
                    for j in natoms:
                        # Look up residue by its number
                        for resi in ob.OBResidueIter(mol):
                            if resi.GetNum() == int(j.split(':')[0]):
                                r = resi
                                break
                        ato = pr.getResAtomByName(r, j.split(':')[1])
                        nnatoms.append(str(ato.GetIdx()))

                    t = "{0[0]:>9s} {0[1]:>9s} {0[2]:>9s} {0[3]:>9s} {1}"
                    r = ""
                    for j in i.split()[4:]:
                        r += "  {}".format(j)
                    r += "  ; {0:>9s} replaced by {1}".format(ra, natoms[nr])
                    marked_impropers.append(t.format(nnatoms, r))
                    # This is marked improper, continue with the next
                    continue
                for b in natoms:
                    for resi in ob.OBResidueIter(mol):
                        # Look up residue by its number
                        if resi.GetNum() == int(b.split(':')[0]):
                            r = resi
                            break
                    ato = pr.getResAtomByName(r, b.split(':')[1])
                    imp.append(ato.GetIdx())
                if len(imp) != 4:
                    s = "Error: Could not find improper {0} in residue {1}"
                    print(s.format(i, res.GetNum()))
                    print(imp)
                si = ""
                for j in i.split()[4:]:
                    si += "     {}".format(j)
                imp.append(si)
                impropers.append(imp)

    impropers = pr.uniquify(impropers)
    impropers.sort()
    marked_impropers = pr.uniquify(marked_impropers)
    marked_impropers.sort()

    return impropers, marked_impropers, hasImpropers


def prefixSugars(mol, n1, n2):
    """Add prefix 'G_' to atomtypes for all atoms in given range.

    """
    s1 = set(range(n1, n2))
    for atom in ob.OBMolAtomIter(mol):
        if (atom.GetIdx() in s1
            and atom.ISYMBL not in ['C', 'H', 'O']):
            atom.ISYMBL = "G_" + atom.ISYMBL


def suffixSugars(mol, n1, n2):
    """Add suffix '_S' to atomtypes for all atoms in given range.

    """
    s1 = set(range(n1, n2))
    for atom in ob.OBMolAtomIter(mol):
        if atom.GetIdx() in s1:
            atom.ISYMBL += "_S"


def writeGromacsTop (mol, itpFN, impTypes={}, sc14={}, tfnTypes=None, defTorsions={}):
    """Write Gromacs topology.

    This function takes an OBMol object and writes the corresponding Gromacs
    topology file (itp).
    """
    nexcl = 3 # minim distance between atoms to have long range interactions
    # Bonded parameter function types
    bfunc = 1 # bond
    pfunc = 1 # 1-4 pair
    afunc = 1 # angle
    dfunc = 9 # dihedral
    ifunc = 2 # improper

    # Default function types for Amber
    fnTypes = OrderedDict([('bfn', 1),
                           ('pfn', 1),
                           ('afn', 1),
                           ('dfn', 9),
                           ('ifn', 4)])

    if tfnTypes is not None:
        for n, t in zip(fnTypes.keys(), tfnTypes):
            if t is not None:
                fnTypes[n] = t
    tbfn, tpfn, tafn, tdfn, tifn = fnTypes.values()

    # Get protein end
    if not mol.HasData('PEND'):
        print("Error: Lost the protein end!")
        sys.exit(4)
    PEND = int(mol.GetData('PEND').GetValue())

    bonds    = pr.getBonds(mol)
    angles   = pr.getAngles(mol)
    torsions = pr.getTorsions(mol)

    pairs14 = pr.getPairs(torsions, angles, bonds)

    # Find impropers in molecule
    imps = pr.findImpropers(mol, impTypes)

    # Get explicit impropers from PREP file if any
    impropers, hasImpropers = pr.getExplicitImps(mol)

    # Add impropers stored in atoms (protein and marked)
    impropers, marked_impropers, hasImpropers = getResidueImps(mol, impropers,
                                                               hasImpropers)

    # Add together found and explicit impropers
    eimps = []
    for i in impropers:
        eimps.append(set(i[:4]))
    for i in marked_impropers:
        eimps.append(set(int(j) for j in i.split()[:4]))
    if imps != []:
        for imp in imps:
            if set(imp[:4]) not in eimps and max(imp[:4]) >= PEND:
                impropers.append(imp)
                impropers[-1][4] = 4
        hasImpropers = True

    # Make protein atom set and sugar atom set
    sp = set(range(1, PEND))
    ss = set(range(PEND, mol.NumAtoms()+1))
    # Mark all lists
    mark_cross(mol, sp, ss, bonds)
    mark_cross(mol, sp, ss, pairs14)
    mark_cross(mol, sp, ss, angles)
    mark_cross(mol, sp, ss, torsions)
    mark_cross(mol, sp, ss, impropers)

    pr.calc14NBParams(mol, pairs14, sc14)

    if defTorsions:
        setDefTorsions(torsions, defTorsions)

    if (options.opls):
        if options.amber:
            options.opls == False
        else:
            print(" Suffixing polymer atomtypes ...")
            suffixSugars(mol, PEND, mol.NumAtoms()+1)
            print(" Done.")

    # Construct the ITP file
    t = []
    t.append("[ moleculetype ]\n")
    t.append("; Name    nrexcl\n")
    props = {"MolNM": "Molecule", "nrexcl": nexcl}
    t.append("{0[MolNM]:s} {0[nrexcl]:d}\n".format(props))

    t.append("\n[ atoms ]\n")
    t.append(";   nr       type  resnr residue  atom  cgnr    charge"
            "     mass\n")
    runningQ  = 0.0
    runningCG = 0
    s1 = "; Residue {0:d} {1:s} q={2:.3f}\n"
    s2 = "{0[aID]:6d}{0[aT]:>11s}{0[resID]:7d}{0[resN]:>7s}" \
         "{0[aN]:>7s}{0[aCG]:7d}{0[aQ]:11.4f}{0[aM]:11.3f}{0[aC]}\n"
    for i in ob.OBResidueIter(mol):
        t.append(s1.format(i.GetNum(), i.GetName(), pr.getResCharge(i)))
        if i.HasData("NOTE"):
            t.append("; {0:s}\n".format(i.GetData("NOTE").GetValue()))
        props = {"resID": i.GetNum(), "resN": i.GetName()}
        for j in ob.OBResidueAtomIter(i):
            runningQ += j.GetPartialCharge()
            props["aID"] = j.GetIdx()
            props["aT"]  = j.ISYMBL
            props["aN"]  = j.IGRAPH
            props["aQ"]  = j.GetPartialCharge()
            props["aM"]  = j.GetAtomicMass()
            props["aC"]  = " ; qtot {0:8.4f}".format(runningQ)
            if j.HasData("CHARGE_GROUP"):
                props["aCG"] = int(j.GetData("CHARGE_GROUP").GetValue())
                runningCG = props["aCG"]
            else:
                runningCG += 1
                props["aCG"] = runningCG
            t.append(s2.format(props))

    # TODO: Print parameters for crosslink?
    # Bonds
    t.append("\n[ bonds ]\n")
    t.append(";   ai     aj   funct   c0   c1\n")
    s = "{0[0]:6d} {0[1]:6d} {1:6d}{0[2]}\n"
    for i in bonds:
        t.append(s.format(i, bfunc))

    # 1-4 pairs
    t.append("\n[ pairs ]\n")
    t.append(";  ai    aj   funct   c0   c1\n")
    s1 = "{0[0]:5d} {0[1]:5d} {1:5d}{2:s}\n"
    s2 = "{0[0]:6d} {0[1]:6d} {1:6d}   {2:6.4f} {3:9.4f} {4:9.4f}" \
         "   {5:1.5e}   {6:1.5e}{7:s}\n"
    for i in pairs14:
        s3 = set(i[:2])
        if not s3.isdisjoint(sp):
            t.append(s1.format(i, pfunc, i[2]))
        else:
            rem = i[2] + " ; Sugar only"
            if len(i) > 3:
                t.append(s2.format(i, 2, i[3], i[4], i[5], i[6], i[7], rem))
            else:
                t.append(s1.format(i, pfunc, rem))

    # Angles
    t.append("\n[ angles ]\n")
    t.append(";    ai    aj     ak   funct   c0   c1\n")
    s = "{0[0]:6d} {0[1]:6d} {0[2]:6d} {1:6d}{0[3]}\n"
    for i in angles:
        t.append(s.format(i, afunc))

    # Proper torsions
    # Print function type from topology for protein and crosslink,
    # and function 9 for sugar
    t.append("\n[ dihedrals ]\n")
    t.append(";   ai     aj     ak     al  funct\n")
    s = "{0[0]:6d} {0[1]:6d} {0[2]:6d} {0[3]:6d} {1:6d} {0[4]}\n"
    for i in torsions:
        s3 = set(i)
        if s3.isdisjoint(sp):
            t.append(s.format(i, dfunc))
        else:
            t.append(s.format(i, tdfn))

    # Improper torsions
    s = "{0[0]:6d} {0[1]:6d} {0[2]:6d} {0[3]:6d} {1:6d}"
    r = "{0[0]:6d} {0[1]:6d} {0[2]:6d} {0[3]:6d}{0[4]}\n"
    if hasImpropers:
        t.append("\n[ dihedrals ]\n")
        t.append(";   ai     aj     ak     al  funct\n")
        for i in impropers:
            if i[4] != 4:
                t.append(r.format(i))
            else:
                if i[4] == 4:
                    t.append(s.format(i[:4], 4) + "\n")
                else:
                    t.append(s.format(i[:4], ifunc) +
                            "    0.000     167.360{0[4]}\n".format(i))
        for i in marked_impropers:
            t.append(" {}\n".format(i))

    t.append("\n")

    with open(itpFN, 'w') as f:
        f.writelines(t)

def charmmGui (inputPDB, inputITP):
    """Write Charmm-gui compatible PDB file.

    This function takes the output PDB and ITP files and writes the corresponding
    Charmm-Gui compatible PDB file.
    """
    atom        = []
    hetatm      = []
    hetatm_old = []
    header      = "TITLE CHARMM-GUI COMPATIBLE\n"
    bondsITP    = []
    boxsize = []
    box = re.compile(r"CRYST")

    with open(options.pathFN  , "r") as f:
        glycamResList = [i.strip("\n") for i in f.readlines()] 

    glycamRes = re.compile(r'\b(?:%s)\b' % '|'.join(glycamResList))


    with open(inputPDB) as file:    
        for line in file:
            if  not re.findall(glycamRes,line):
                if line.startswith("ATOM"):
                    atom.append(line)
                if line.startswith("HETATM"):
                    hetatm_old.append(line)

            elif re.findall(glycamRes,line):
                    hetatm.append(line)
            if re.match(box,line):
                    boxsize.append(line)

    with open(inputITP, "r") as f:
        bondsITP.append(getITPBonds(f))


    bonds= []
    for i in range(len(bondsITP[0])):
        if (bondsITP[0][i][0] >= len(atom)+1):
                bonds.append(bondsITP[0][i])    

        if (bondsITP[0][i][1] >= len(atom)+1):
                bonds.append(bondsITP[0][i]) 
    bonds = pr.uniquify(bonds)
    
    # Replace ATOM with HETATM and viceversa:
    hetatm = [w.replace("ATOM  " , "HETATM") for w in hetatm]
    atom_new    = [w.replace("HETATM" , "ATOM  ") for w in hetatm_old]

        
    output = str(inputPDB).strip(".pdb") + "_charmm-gui.pdb"
    
    with open(output,"w") as newfile:    
        newfile.write(header)
        newfile.write("".join(boxsize))
        newfile.write("".join(atom))
        newfile.write("".join(atom_new))
        newfile.write("TER\n")
        newfile.write("".join(hetatm))
        for i in bonds:
            newfile.write ("%6s%5s%6s" % ("CONECT", str(i[0]), str(i[1]) + "\n" ))
        newfile.write("END")


if __name__ == "__main__":
    options, args = optP()

    print("Reading forcefield ...")
    if  not os.path.isfile(options.ffFN):
        print("Error: File {} not found!".format(options.ffFN))
        sys.exit(1)
    sc14 = pr.get14NBParams(options.ffFN)
    imps = pr.getImproperTypes(options.ffFN)
    print("Done.")

    print("Reading topology ...")
    with open(options.itpFNin, 'r') as f:
        # Get atom types, impropers and function types from topology
        atypes, resnumbers, resname, anames, chgroups, charges = getAtomTypes(f)
        impropers, improper_atoms = getImpropers(f)
        fnTypes = getFunctionTypes(f)
        bonds = getITPBonds(f)
        defTorsions = getDefTorsions(f)
    print("Done.")
    # Recode defined torsions using atom names and residue numbers
    s = "{}:{}"
    defTorsions2 = {}
    for i in defTorsions.keys():
        an = s.format(resnumbers[i[0]], anames[i[0]])
        bn = s.format(resnumbers[i[1]], anames[i[1]])
        cn = s.format(resnumbers[i[2]], anames[i[2]])
        dn = s.format(resnumbers[i[3]], anames[i[3]])
        defTorsions2[(an,bn,cn,dn)] = defTorsions[i]
    defTorsions.update(defTorsions2)
    # Recode impropers using atom names
    # Take into account residue number for atom
    impstr = " {0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:7d}    {5}"
    ifa = {i: '' for i in improper_atoms}
    for i in impropers:
        an = s.format(resnumbers[i[0]], anames[i[0]])
        bn = s.format(resnumbers[i[1]], anames[i[1]])
        cn = s.format(resnumbers[i[2]], anames[i[2]])
        dn = s.format(resnumbers[i[3]], anames[i[3]])
        if len(i) == 5:
            impS = impstr.format(an, bn, cn, dn, i[4], '')
        else:
            impS = impstr.format(an, bn, cn, dn, i[4], i[5][0])
        ifa[i[0]] += "{}-".format(impS)
        ifa[i[1]] += "{}-".format(impS)
        ifa[i[2]] += "{}-".format(impS)
        ifa[i[3]] += "{}-".format(impS)

    print("Reading PREP file(s) ...")
    residues = {}
    for i in options.prepFNs:
        print(i)
        with open(i, 'r') as f:
            ires = pr.readPrepFile(f, options.no_warn)
        for k, v in ires.items():
            residues[k] = v

    print("Got", len(residues), "residues in library")

    print("Reading protein structure ...")
    # Read in provided structure file and parse molecules in it
    fformat = os.path.splitext(options.protFN)[1][1:]
    if not options.noBonds:
        pybmol = next(pybel.readfile(fformat, options.protFN))

        mol = pybmol.OBMol
        mol.SetAutomaticPartialCharge(False)
    else:
        mol = ob.OBMol()
        mol.SetAutomaticPartialCharge(False)

        conv = ob.OBConversion()
        conv.SetInFormat(fformat)
        # Add option to disable bond perception on reading
        # Relies on connection table for bonding information
        conv.AddOption('b', conv.INOPTIONS)
        conv.ReadFile(mol, options.protFN)

        # Clear all bonds if any (e.g., from CONECT records in PDB file)
        old_bonds = [bond for bond in ob.OBMolBondIter(mol)]
        for bond in old_bonds:
            mol.DeleteBond(bond)

        # Add bonds from ITP file
        for a, b in bonds:
            pr.addBond(mol, mol.GetAtom(a), mol.GetAtom(b))

        # Still need to perceive bond orders for SMARTS to work as expected
        mol.PerceiveBondOrders()

        pybmol = pybel.Molecule(mol)
    print("Got {} atoms".format(len(pybmol.atoms)))

    print("Molecule has {} bonds".format(mol.NumBonds()))
    print("Molecule has {} residues".format(mol.NumResidues()))

    print("Marking residues and atoms ...")
    for res in ob.OBResidueIter(mol):
        resAts = set()
        resNr  = res.GetNum()
        # OB defaults to PDB 3 letter residue names
        # while Gromacs can use 4 letter residue names
        # set it explicitly to preserve naming in topology
        res.SetName(resname[resNr])
        # Add chain IDs to residues
        pr.addPairData(res, 'TRUNK_ID',  "P.{}".format(resNr-1))
        pr.addPairData(res, 'BRANCH_ID', "P.{}".format(resNr))
        # Add "IGRAPH" data for protein atoms
        for atom in ob.OBResidueAtomIter(res):
            res.SetAtomID(atom, anames[atom.GetIdx()])
            atom.IGRAPH = anames[atom.GetIdx()]
            atom.ISYMBL = atypes[atom.GetIdx()]
            atom.CHRG   = charges[atom.GetIdx()]
            pr.addPairData(atom, 'CHARGE_GROUP', chgroups[atom.GetIdx()])
            atom.SetPartialCharge(float(charges[atom.GetIdx()]))
            resAts.add(atom.GetIdx())
        common = resAts.intersection(improper_atoms)
        if common:
            impS = ""
            for a in common:
                impS += ifa[a]
            pr.addPairData(res, 'IMPROPERS', impS.rstrip('-'))
    print("Done.")

    print("Reading sequence file(s) ...")
    sequence = []
    nreplace = 0
    reXr = re.compile(r'^\w+=P[.][0-9]+')
    with open(options.seqFN, 'r') as f:
        lines = f.readlines()
        for i in lines:
            if reXr.match(i) is not None:
                nreplace += 1
            sequence.append(i.rstrip('\n'))
    print(" Found {} replaced atoms".format(nreplace))
    # Protein end = size - number of replaced atoms + 1
    pr.addPairData(mol, 'PEND', str(mol.NumAtoms() - nreplace + 1))
    print("Done.")

    print("Glycosylating protein ...")
    omol = buildSeq(residues, mol, sequence)
    print("Done.")

    print("Writing Gromacs topology ...")
    writeGromacsTop(omol, options.itpFNout, imps, sc14, fnTypes, defTorsions)
    print("Done.")

    print("Writing structure file ...")
    opybmol = pybel.Molecule(omol)
    offormat = os.path.splitext(options.postFN)[1][1:]
    opybmol.write(offormat, options.postFN, overwrite=True)
    print("Done.")

    
    if (options.pathFN  ):
        print ("Writing Charmm-Gui compatible PDB file")
        if offormat != "pdb":
            print ("\nWARNING: your output File has to be in PDB format!\n"
                   "The Charmm-Gui compatible PDB file has not been created") 
        else :
            charmmGui(options.postFN,options.itpFNout)
            print ("\nNOTE : The Charmm-Gui won't read structures with steric clashes.\n" \
               "Please, consider to minimize your output structure first and then run\n" \
               "doglycans2charmmgui.py separately.")
            print ("Done.")
