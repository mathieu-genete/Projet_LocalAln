#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:51:56 2024

@author: mathieu
"""

class IUPAC:
    
    @staticmethod
    def Iupac_code():
        code={'A':{'A'},'C':{'C'},'G':{'G'},'T':{'T'},'R':{'A','G'},
              'Y':{'C','T'},'S':{'G','C'},'W':{'A','T'},'K':{'G','T'},
              'M':{'A','C'},'B':{'C','G','T'},'D':{'A','G','T'},
              'H':{'A','C','T'},'V':{'A','C','G'},'N':{'A','C','G','T'}}
        return code