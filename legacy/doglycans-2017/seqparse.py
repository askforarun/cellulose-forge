#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Version 0.01
Author: See the AUTHORS file for a list of contributors to doGlycans.py
License GPLv3

Parser for sequence files.
"""

import sys
import ply.yacc as yacc

from seqtokens import tokens


def p_branch9(p):
    'branch : base connection segments TIP'
    p[0] = p[1] + [p[2]] + p[3] + [p[4]]

def p_branch8(p):
    'branch : base connection TIP'
    p[0] = p[1] + [p[2], p[3]]

def p_branch7(p):
    'branch : base connection NAME tag'
    p[0] = p[1] + [p[2], p[3], tuple([None, None] + [i for i in p[4]])]

def p_branch6(p):
    'branch : root attach bond COMMA conf'
    p[0] = p[1] + [p[2]] + [tuple([p[3], p[5]])]

def p_branch5(p):
    'branch : root attach NAME COMMA TIP'
    p[0] = p[1] + [p[2]] + [(tuple([p[3], p[5]]), None)]

def p_branch4(p):
    'branch : root segments TIP'
    p[0] = p[1] + [None, None] + p[2] + [p[3]]

def p_branch3(p):
    'branch : root TIP'
    p[0] = p[1] + [None, None, p[2]]

def p_branch2(p):
    'branch : root segments NAME tag'
    p[0] = p[1] + [None, None] + p[2] \
                + [p[3], tuple([None, None] + [i for i in p[4]])]

def p_branch1(p):
    'branch : root bond COMMA conf'
    p[0] = p[1] + [None] + [tuple([p[2]] + [p[4]])]

def p_base2(p):
    'base : root attach'
    p[0] = p[1] + [p[2]]

def p_base1(p):
    'base : root'
    p[0] = p[1] + [None]

def p_root(p):
    'root : NAME EQUAL'
    p[0] = [p[1]]

def p_attach2(p):
    'attach : ATTACH REPLACE NAME COLON'
    p[0] = (p[1], p[3])

def p_attach1(p):
    'attach : ATTACH COLON'
    p[0] = (p[1], None)

def p_block2(p):
    'segment : block'
    p[0] = p[1]

def p_block1(p):
    'block : LCBR segments RCBR NUMBER'
    p[0] = p[2] * p[4]

def p_segments2(p):
    'segments : segments segment'
    p[0] = p[1] + p[2]

def p_segments1(p):
    'segments : segment'
    p[0] = p[1]

def p_segment2(p):
    'segment : NAME tag connection'
    # Put all the labels as strings after angles
    p[0] = [p[1], tuple(list(p[3]) + p[2])]

def p_segment1(p):
    'segment : NAME connection'
    p[0] = [p[1], p[2]]

def p_connection3(p):
    'connection : LPAREN bond COMMA tag RPAREN'
    # Only one angle label per connection can be used
    p[0] = tuple([p[2]] + p[4][:1])

def p_connection2(p):
    'connection : LPAREN bond COMMA conf RPAREN'
    p[0] = (p[2], p[4])

def p_connection1(p):
    'connection : LPAREN bond RPAREN'
    p[0] = (p[2], None)

def p_tag(p):
    'tag : LABR labels RABR'
    p[0] = p[2]

def p_labels2(p):
    'labels : labels LABEL'
    p[0] = p[1] + [p[2]]

def p_labels1(p):
    'labels : LABEL'
    p[0] = [p[1]]

def p_bond(p):
    'bond : NAME COMMA NAME'
    p[0] = (p[1], p[3])

def p_conf(p):
    'conf : LSBR angles RSBR'
    p[0] = p[2]

def p_angles2(p):
    'angles : angles angle'
    p[0] = p[1] + [p[2]]

def p_angles1(p):
    'angles : angle'
    p[0] = [p[1]]

def p_angle(p):
    'angle : NAME NAME NAME NAME ANGLE'
    p[0] = (p[1], p[2], p[3], p[4], p[5])

def p_error(p):
    print("Syntax error in input!")
    sys.exit(1)


parser = yacc.yacc()


if __name__ == '__main__':
    s1 = "HM=P.150/HD22:-(ND2,C1)-{4YB-(O4,C1)-4YB<C3,C4>-(O4,C1,<Con,C5>)}3" \
         "-VMB-(O6,C1,[H5 C5 C6 O6 0.0 C5 C6 O6 C1 180.0 C6 O6 C1 H1 60.0])-0MA"
    s2 = "branch1=HM.3:-(O3,C1,[C3 O3 C1 H1   0.0])-0MA"
    s3 = "Mol=ROH-(O1,C1)-{4YB-(O4,C1)-4YB-(O4,C1,[H4 C4 O4 C1 -15.0])}2-0MA"
    s4 = "0XU=ROH-(O1,C1)-0XU"
    s5 = "branch2=HM.3:-(O6,C1,<Con,C3>)-0MA"
    l1 = "Con=O6,C1,[H5 C5 C6 O6 0.0 C5 C6 O6 C1 180.0 C6 O6 C1 H1 60.0]"
    a1 = "HM1=P.150/HD22:-(ND2,C1)-4YB-(O4,C1)-4YB-" \
         "(O4,C1,<HM1a>)-VMB<HM1b1>-(O6,C1,<HM1b>)-0MA"
    b1 = "Bon=P.96:ND1,FE"
    b2 = "Bon=P.96:ND1,FE,[C3 O3 C1 H1   0.0 C3 O3 C1 H1   0.0]"
    r1 = "0XU=0XU"
    print(parser.parse(s1, debug=1))
