#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 15:48:34 2024

@author: mathieu
"""
from .Alphabet import Alphabet

class Seq:
    
    def __init__(self,id,seq,description="",alphabet=Alphabet.IUPAC()):
        self.__alphabet=alphabet
        self.__id=id
        if description=="":
            self.__description=id
        else:
            self.__description=description
        self.__seq=seq.upper()
        self.__check_dna()
        
    #===================
    #Méthodes magiques
    #=================== 
    
    def __str__(self):
        return self.__seq
    
    #===================
    #Getters Setters
    #===================
    
    @property
    def alphabet(self):
        return self.__alphabet
    
    @property
    def seq(self):
        return self.__seq
    
    @property
    def description(self):
        return self.__description
    
    @property
    def id(self):
        return self.__id
        
    #===================
    #Méthodes privées
    #===================
    def __check_dna(self):
        for b in self.__seq:
            if b not in self.__alphabet:
                raise Exception("Error alphabet in sequence id: {}".format(self.__id))