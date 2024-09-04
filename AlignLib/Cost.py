#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 21 19:55:25 2024

@author: mathieu
"""
from .Alphabet import Alphabet
from .IUPAC import IUPAC

class Cost:
    
    REQ_COSTS_PARAM={'identity','substitution','indel','gap_open','gap_ext','alternative_match'}
    
    def __init__(self,scores: dict,alphabet=Alphabet.IUPAC()):
                
        if not Cost.REQ_COSTS_PARAM.issubset(set(scores.keys())):
            req_param=list(Cost.REQ_COSTS_PARAM.difference(set(scores.keys())))
            raise Exception("Il manque des paramètres dans le fichier de configuration: {}".format(" - ".join(req_param)))
        
        self.__alphabet=alphabet
        self.__identity = scores['identity']
        self.__substitution = scores['substitution']
        self.__indel = scores['indel']
        self.__gap_open = scores['gap_open']
        self.__gap_ext = scores['gap_ext']
        self.__alternative_match = scores['alternative_match']
        self.__infini = float('inf')
                
    def intersect_IUPAC(self,b1,b2):
        code=IUPAC.Iupac_code()
        return len(code[b1].intersection(code[b2]))
        
    def __compare(self,b1,b2):
        if b1 not in self.__alphabet or b2 not in self.__alphabet:
            print(b1,b2)
            raise Exception("Opération invalide pour ces caractères.")
            
        if self.__alphabet == Alphabet.dna():
            return self.__costACGT(b1,b2)
        elif self.__alphabet == Alphabet.IUPAC():
            if b1 in Alphabet.dna() and b2 in Alphabet.dna():
                return self.__costACGT(b1,b2)
            else:
                intersect=self.intersect_IUPAC(b1,b2)
                if intersect==0:
                    return self.__substitution
                else:
                    return self.__alternative_match
                
    def __costACGT(self,b1,b2):
        if b1==b2:
            return self.__identity
        elif b1!=b2:
            return self.__substitution
           
    #===================
    #Méthodes magiques
    #=================== 
    
    def __call__(self,b1,b2):
        return self.__compare(b1,b2)
    
    def __str__(self):
        out_info ="Alphabet: {}\n".format(self.__alphabet)
        out_info+="Costs:\n"
        out_info+="  Indentity: {}\n".format(self.__identity)
        out_info+="  Substitution: {}\n".format(self.__substitution)
        out_info+="  Indel: {}\n".format(self.__indel)
        out_info+="  Gap_penalty: {}\n".format(self.__gap_open)
        out_info+="  Extend_penalty: {}\n".format(self.__gap_ext)
        return out_info
    
    #===================
    #Getters Setters
    #===================
    
    @property
    def indel(self):
        return self.__indel
    
    @property
    def identity(self):
        return self.__identity
    
    @property
    def substitution(self):
        return self.__substitution
    
    @property
    def gap_open(self):
        return self.__gap_open
    
    @property
    def gap_ext(self):
        return self.__gap_ext
    
    @property
    def infini(self):
        return self.__infini