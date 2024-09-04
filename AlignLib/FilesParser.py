#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:45:10 2024

@author: mathieu
"""

class FilesParser:
    
    @staticmethod
    def ParseConfig(file):
        config={}
        with open(file,'r') as inconfig:
            for line in inconfig:
                line=line.strip()
                tmp=line.split("#")[0]
                if len(tmp)>0:
                    key_value=tmp.split(":")
                    config[key_value[0].strip()]=FilesParser.__check_value(key_value[1])
        return config
    
    def __check_value(value):
        if FilesParser.__is_int(value):
            return int(value)
        elif FilesParser.__is_float(value):
            return float(value)
        else:
            return value
        
    def __is_int(n):
        try:
            float_n = float(n)
            int_n = int(float_n)
        except ValueError:
            return False
        else:
            return float_n == int_n

    def __is_float(n):
        try:
            float_n = float(n)
        except ValueError:
            return False
        else:
            return True
