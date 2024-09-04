#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Programme d'alignement local de séquences
Usage:
======
python localAln.py -h
pour afficher l'aide

@author: Mathieu Genete
"""
__authors__ = ("Mathieu GENETE")
__contact__ = ("mathieu.genete@univ-lille.fr")
__copyright__ = "copyleft"
__date__ = "2024/01"
__version__= "1.0.0"

import argparse
import sys
from AlignLib import FastaSeq, FilesParser, Cost, Alignment

def main():
    
    #Définition des paramètres du programme avec argparse
    description=""" Programme d'alignement local de séquences """
    parser=argparse.ArgumentParser(prog="localAln.py",description=description)
    parser.add_argument("-i","--infasta", help="fasta contenant les séquences à analyser", required = True)
    parser.add_argument("-c","--config", help="fichier de configuration defaut='datas/score.txt'", default="datas/score.txt")
    parser.add_argument("-o","--outstats", help="exporter les statistiques d'alignement")
    parser.add_argument("-f","--outfasta", help="exporter l'alignement au format fasta")
    parser.add_argument("-l","--linesize", help="nombre de bases par ligne pour l'alignement",type=int,default=50)
    parser.add_argument("-n","--rndseq", help="nombre de séquences aléatoires à générer pour le calcul de la significativité (n=100 par défaut)",type=int,default=100)
    parser.add_argument("-g","--gapaffine", help="Utiliser le système de score affine pour les indels",action='store_const',const=True,default=False)
    parser.add_argument("-s","--sinificiance", help="Calcul de la significativité du score d'alignement",action='store_const',const=True,default=False)
    parser.add_argument("-M","--showmatrix", help="Afficher la matrice d'alignement",action='store_const',const=True,default=False)

    args = parser.parse_args(sys.argv[1:])
    
    fastafile=args.infasta
    
    fasta = FastaSeq(fastafile)
    if len(fasta.seqlist)>=2:
        U=fasta.seqlist[0]
        V=fasta.seqlist[1]
    
    scores = FilesParser.ParseConfig(args.config)    
    cost = Cost(scores)
    
    aln1 = Alignment(U,V,cost,GapAffine=args.gapaffine,lineSize=args.linesize)
    
    aln1.Display_Align()
    
    if args.showmatrix:
        aln1.print_matrix()
    
    if args.sinificiance:
        aln1.Significance(args.rndseq)
        aln1.DisplayDist()
        
    if args.outstats:
        aln1.Export_Alignment(args.outstats,outtype="stats")
    
    if args.outfasta:
        aln1.Export_Alignment(args.outfasta,outtype="fasta")
     
    
if __name__=="__main__":
    main()