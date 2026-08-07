#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Version 0.01
Author: See the AUTHORS file for a list of contributors to doGlycans.py
License GPLv3

Module with lexing rules for sequence files.
"""

import ply.lex as lex


states = (
    ('angle',  'exclusive'),
    ('square', 'exclusive'),
)

tokens = (
    'ANGLE',
    'ATTACH',
    'COLON',
    'COMMA',
    'EQUAL',
    'NAME',
    'NUMBER',
    'REPLACE',
    'LABEL',
    'LPAREN',
    'RPAREN',
    'LSBR',
    'RSBR',
    'LCBR',
    'RCBR',
    'LABR',
    'RABR',
    'TIP',
)

t_COLON   = r':'
t_COMMA   = r','
t_EQUAL   = r'='

t_REPLACE = r'/'
t_LPAREN  = r'\('
t_RPAREN  = r'\)'
t_LCBR    = r'\{'
t_RCBR    = r'\}'

def t_ATTACH(t):
    r'[a-zA-Z][a-zA-Z0-9_]*[.]\d+'
    return t

def t_TIP(t):
    r'[\d]?[a-zA-Z][a-zA-Z0-9]*\s*$'
    t.value.rstrip()
    return t

def t_NAME(t):
    r'[\d]?[a-zA-Z][a-zA-Z0-9_]*'
    return t

def t_NUMBER(t):
    r'\d+'
    t.value = int(t.value)
    return t

def t_LABR(t):
    r'\<'
    t.lexer.push_state('angle')
    return t

t_angle_LABEL = r'[a-zA-Z][.a-zA-Z0-9_]*'

def t_angle_RABR(t):
    r'\>'
    t.lexer.pop_state()
    return t

def t_LSBR(t):
    r'\['
    t.lexer.push_state('square')
    return t

t_square_NAME = r'[a-zA-Z][a-zA-Z0-9]*'

def t_square_ANGLE(t):
    r'[+-]?\d+[.]?\d*'
    t.value = float(t.value)
    return t

def t_square_RSBR(t):
    r'\]'
    t.lexer.pop_state()
    return t

t_ignore  = ' \t-'
t_angle_ignore = ' \t,'
t_square_ignore = ' \t'

def t_ANY_error(t):
    print("Illegal character {}".format(t.value[0]))
    t.lexer.skip(1)


if __name__ == '__main__':
    lexer = lex.lex(debug=1)
else:
    lexer = lex.lex()
