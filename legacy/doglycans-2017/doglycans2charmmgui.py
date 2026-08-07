#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Version 1
Authors: See the AUTHORS file for a list of contributors to doGlycans.py
License GPLv3

Write a Charmm-Gui compatible PDB file for Glycoproteins strutures. 


Input:
    - Glycosylated protein structure (PDB)
    - Glycosylated protein topology (for reading bonding)
    - Glycam file with sugar residue names (e.g $DGPATH/dat/GLYCAM_06h_resname.txt)

Output:
    - Charmm-Gui compatible PDB file

Note: 
    - The script modifies the PDB file to be read by Charmm-Gui web server.
      It adds the connection records from the topology file and lists the 
      non sugar residues as "ATOM" and sugar residues as "HETATM".
      Please, consider to add the chain ID which can be added using GROMACS
      tool (e.g gmx editconf -f ProtG.pdb -o ProtG.pdb -label A).
      For more information on Charmm-Gui compatible files, please follow the below link:
      http://www.charmm-gui.org/?doc=input/glycan
"""

import re
import sys
import os.path
import numpy as np
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

    usage="[python3] %prog -c ProtG.pdb -i ProtG.itp --charmm" \
          " $DGPATH/dat/GLYCAM_06h_resname.txt"

    description=""
    version="\n%prog Version 0.04\n\nRequires Python 3, Open Babel Python" \
            " bindings, numpy and prepReader.py."

    optParser = OptionParser(usage=usage,
                             version=version,
                             description=description)

    optParser.set_defaults(ffFN=None,
                        )

    
    optParser.add_option('-i', type='str',
                         dest='itpFN',
                         help="Gromacs topology file (in) [default: %default]")
    optParser.add_option('-c', type='str',
                         dest='protFN',
                         help="Protein structure file [default: %default]")
    optParser.add_option('--charmm', type='str',
                         dest='pathFN',
                         help="Print PDB for charmm-gui submission.\n" \
                         "The path to a file with sugar residue name list is needed"
                          " [default: %default]")
    
    options, args = optParser.parse_args()



    if options.protFN is None:
        s = "Error: PDB structure file required."
        print(s)
        sys.exit(2)

    if options.itpFN is None:
        s = "Error: Input .itp file required."
        print(s)
        sys.exit(2)

    if options.pathFN is None:
        s = "Error: Glycam residue name list required."
        print(s)
        sys.exit(2)
    return options, args

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

    with open(options.pathFN, "r") as f:
        glycamResList = [i.strip("\n") for i in f.readlines()] 

    glycamRes = re.compile(r'\b(?:%s)\b' % '|'.join(glycamResList))
    box = re.compile(r"CRYST")

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

    print ("Writing Charmm-Gui compatible PDB file")
    if os.path.splitext(options.protFN)[-1] != ".pdb":
        print ("\nWARNING: your output File has to be in PDB format!\n"
               "The Charmm-Gui compatible PDB file has not been created") 
    else :
        charmmGui(options.protFN,options.itpFN)
        print ("\nNOTE : The Charmm-Gui won't read structures with steric clashes.\n" \
               "Please, consider to minimize your output structure first and then run\n" \
               "doglycans2charmmgui.py separately.")
        print ("Done.")