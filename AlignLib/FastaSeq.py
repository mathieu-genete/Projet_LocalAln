#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:12:26 2024

@author: mathieu
"""
from .Alphabet import Alphabet
from .Seq import Seq

import os

class FastaSeq:
    def __init__(self,infasta,alphabet=Alphabet.IUPAC()):
        self.__Alphabet=alphabet
        self.__seqlist=[]
        self.__ParseFasta(infasta)        
    
    #===================
    #Méthodes privées
    #===================
    def __ParseFasta(self,infasta):
        try:
            if not os.path.exists(infasta):
                raise Exception(f"Le fichier '{infasta}' n'existe pas")
                
            with open(infasta,'r') as fasta:
                seq={'id':"",'description':"",'seq':""}
                for line in fasta:
                    line=line.strip()
                    if line.startswith(">"):
                        if seq['id']!="":
                            self.__seqlist.append(Seq(id=seq['id'],seq=seq['seq'],description=seq['description'],alphabet=self.__Alphabet))
                            seq={'id':"",'description':"",'seq':""}
                            
                        seq['id']=line[1:].split()[0]
                        seq['description']=line[1:]
                    else:
                        seq['seq']+=line
                self.__seqlist.append(Seq(id=seq['id'],seq=seq['seq'],description=seq['description'],alphabet=self.__Alphabet))
        except:
            raise Exception(f"Erreur avec le fichier '{infasta}'")
            
    #===================
    #Getters Setters
    #===================
    
    @property
    def seqlist(self):
        return self.__seqlist
    
#    @seq.setter
#    def seq(self,value):
#        self.__seq=value
    
    @property
    def id(self):
        return self.__id
        