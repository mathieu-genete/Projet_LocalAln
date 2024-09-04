#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 08:51:33 2024

@author: mathieu
"""
import sys

class Output_screen:
    
    @staticmethod
    def sdterr_print(value=""):
        sys.stderr.write("{}\n".format(value))
        
    @staticmethod
    def sdtout_print(value):
        sys.stdout.write("{}\n".format(value))