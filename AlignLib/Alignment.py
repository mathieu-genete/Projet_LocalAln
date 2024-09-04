#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:20:12 2024

@author: mathieu
"""
import numpy as np
import pandas as pd
import random
from datetime import datetime
from .Output_screen import Output_screen
from AlignLib import __title__, __version__
from .Cost import Cost
from .Seq import Seq

class Alignment:
    
    def __init__(self,SeqU: Seq,SeqV: Seq,cost: Cost,GapAffine=False,mode="",lineSize=80):
        self.__SeqU = SeqU
        self.__SeqV = SeqV
        self.__cost = cost
        self.__sorted_scores = None
        self.__GapAffine = GapAffine
        self.__mode = mode
        self.__lineSize = lineSize
        
        self.__aln_stats = {'identity':0,'similarity':0,'gaps':0,'length':0}
        
        if self.__GapAffine:
            self.__matrix,self.__matrix_I,self.__matrix_D,self.__max_index_T,self.__aln_score = self.__DPmatrix_Gap()
            self.__outAln,self.__min_index = self.__align_traceback_gap()
        else:
            self.__matrix,self.__max_index,self.__aln_score = self.__DPmatrix()
            if self.__mode=="rec":
                i=self.__max_index[0]
                j=self.__max_index[1]
                self.__outAln,self.__min_index = self.__align_traceback_rec(i,j,[])
            else:
                self.__outAln,self.__min_index = self.__align_traceback()
    
        self.__SeqU_Aligned,self.__SeqUV_Matchs,self.__SeqV_Aligned =self.__StoreAlignment(self.__outAln)
        
        self.__outStatsAln = self.__print_Aln(self.__outAln,self.__min_index,self.__lineSize)
    
    #===================
    #Méthodes publiques
    #===================
    
    def print_matrix(self):
        Output_screen.sdterr_print("================== Matrice d'alignement ==========================")
        U=" "+ self.__SeqU.seq
        V=" "+ self.__SeqV.seq
        print_matrix=np.array(self.__matrix)
        Output_screen.sdterr_print(pd.DataFrame(print_matrix.transpose(), index = [*V], columns = [*U]))
        Output_screen.sdterr_print("==================================================================\n")
        
    def Display_Align(self,outtype="stats"):
        if self.__GapAffine:
            Output_screen.sdterr_print("-- use GapAffine (stack traceback) --\n")
        else:
            if self.__mode=="rec":
                Output_screen.sdterr_print("-- use recursive traceback --\n")
            else:
               Output_screen.sdterr_print("-- use loop traceback --\n")
        if outtype=="stats":
            Output_screen.sdtout_print(self.__outStatsAln.format(report_file="stdout"))
        elif outtype=="fasta":
            Output_screen.sdtout_print(self.__AlnToFasta(self.__outAln))

        
    def Significance(self,n=100):
        SeqU_list=list(self.__SeqU.seq)
        scores=[]
        for _ in range(n):
            random.shuffle(SeqU_list)
            rnd_SeqU="".join(SeqU_list)
            if not self.__GapAffine:
                scores.append(self.__DPmatrix_maxscore(rnd_SeqU,self.__SeqV.seq))
            else:
                T,I,D,max_index,max_score = self.__DPmatrix_Gap(rnd_SeqU,self.__SeqV.seq)
                scores.append(int(max_score))
            
        sorted_scores=sorted(scores,reverse=True)
        for sig in range(0,len(sorted_scores)):
            if sorted_scores[sig]<self.__aln_score:
                break
        sig_pct=(float(sig)/float(n))*100
        Output_screen.sdterr_print("Significiance of the score: {:3.2f}% ({} sequences)\n".format(sig_pct,n))
        
        self.__sorted_scores=sorted_scores

    def DisplayDist(self,barsize=50):
        if self.__sorted_scores is not None:
            distribution=self.__GetDist(self.__sorted_scores)
            Output_screen.sdterr_print("\t{:4s}\t{}".format("Score","  #"))
            maxcount=max(distribution,key=lambda x:x[1])[1]
            for score,count in distribution:
                aln_star=""
                if score==self.__aln_score:
                    aln_star="   *"
                bars="="*int((count*barsize)/maxcount)
                Output_screen.sdterr_print("{:2s}\t{:5d}\t{:3d}:{}".format(aln_star,score,count,bars))
        else:
            Output_screen.sdterr_print("Launch Significance() method before to display distribution")
    
    def Export_Alignment(self,outfile,outtype="fasta"):
        Output_screen.sdterr_print("Export Alignment as '{}' format in file : '{}'".format(outtype,outfile))
        if outtype=="fasta":
            with open(outfile,'w') as outfasta:
                fastaTXT=self.__AlnToFasta(self.__outAln)
                outfasta.write(fastaTXT)
                
        elif outtype=="stats":
            with open(outfile,'w') as outstats:
                outstats.write(self.__outStatsAln.format(report_file=outfile))
    
    
    #===================
    #Méthodes privées
    #===================
    
    def __StoreAlignment(self,outAln):
        SeqU_Aln=[]
        SeqUV_Match=[]
        SeqV_Aln=[]
        for b1,al,b2 in outAln:
            SeqU_Aln.append(b1)
            SeqUV_Match.append(al)
            SeqV_Aln.append(b2)
            
        return "".join(SeqU_Aln),"".join(SeqUV_Match),"".join(SeqV_Aln)
            
    def __AlnToFasta(self,outAln):
            
        SeqU_print=self.__MultipleLinePrint(self.__SeqU_Aligned,self.__lineSize)
        SeqV_print=self.__MultipleLinePrint(self.__SeqV_Aligned,self.__lineSize)
        
        outfasta=">{seq1ID}\n{seq1Seq}\n>{seq2ID}\n{seq2Seq}".format(seq1ID=self.__SeqU.description,seq1Seq=SeqU_print,seq2ID=self.__SeqV.description,seq2Seq=SeqV_print)
        
        return outfasta
    
    def __MultipleLinePrint(self,Seqlist,lineSize):
        out_print=""
        for x,v in enumerate(Seqlist):
            if (x+1)%lineSize==0:
                out_print+='\n'
            else:
                out_print+=v
        return out_print
    
    def __GetDist(self,scores):
        dist=[]
        maxScore=max(scores)
        for i in range(0,maxScore+1):
            count=scores.count(i)
            dist.append((i,count))
        return dist

#    def __DPmatrix_maxscore_NP(self,U,V):
#        U=" "+ U
#        V=" "+ V
#        max_score=0
#        M=np.array([[0]*len(V) for _ in range(2)])
#        for i,j in itertools.product(range(1,len(U)),range(1,len(V))):                     
#            c1=M[0,j-1]+self.__cost(U[i],V[j])
#            c2=M[0,j]+self.__cost.indel
#            c3=M[1,j-1]+self.__cost.indel
#            M[1,j]=max(c1,c2,c3,0)
#            if M[1,j]>=max_score:
#                max_score=M[1,j]
#            if j==len(V)-1:
#                M = np.delete(M,0,0)
#                M = np.vstack([M,[0]*len(V)])
#        return max_score
    
    def __DPmatrix_maxscore(self,U,V):
        U=" "+ U
        V=" "+ V
        max_score=0
        M=[[0]*len(V) for _ in range(2)]
        for i in range(1,len(U)):
            for j in range(1,len(V)):                      
                    c1=M[0][j-1]+self.__cost(U[i],V[j])
                    c2=M[0][j]+self.__cost.indel
                    c3=M[1][j-1]+self.__cost.indel
                    M[1][j]=max(c1,c2,c3,0)
                    if M[1][j]>=max_score:
                        max_score=M[1][j]
            M=M[1:]+[[0]*len(V)]
        return max_score

    def __DPmatrix_NP(self):
        U=" "+ self.__SeqU.seq
        V=" "+ self.__SeqV.seq
        M=np.zeros((len(U),len(V)))
        max_index=(0,0)
        max_score=0
        for i in range(1,len(U)):
            for j in range(1,len(V)):   
                    c1=M[i-1,j-1]+self.__cost(U[i],V[j])
                    c2=M[i-1,j]+self.__cost.indel
                    c3=M[i,j-1]+self.__cost.indel
                    M[i,j]=max(c1,c2,c3,0)
                    if M[i,j]>=max_score:
                        max_score=M[i,j]
                        max_index=(i,j)
        return M,max_index,max_score
    
    def __DPmatrix(self):
        U=" "+ self.__SeqU.seq
        V=" "+ self.__SeqV.seq
        M=[[0]*(len(V)) for _ in range(len(U))]
        max_index=(0,0)
        max_score=0
        for i in range(1,len(U)):
            for j in range(1,len(V)):
                c1=M[i-1][j-1]+self.__cost(U[i],V[j])
                c2=M[i-1][j]+self.__cost.indel
                c3=M[i][j-1]+self.__cost.indel
                M[i][j]=max(c1,c2,c3,0)
                if M[i][j]>max_score:
                    max_score=M[i][j]
                    max_index=(i,j)
        return M,max_index,max_score
    
    def __DPmatrix_Gap(self,U="",V=""):
        if U=="":
            U=" "+ self.__SeqU.seq
        else:
            U=" "+ U
            
        if V=="":
            V=" "+ self.__SeqV.seq
        else:
            V=" "+ V
                    
        #initialisation matrices
        T=[[0]*len(V) for _ in range(len(U))]
        I=[[-self.__cost.infini]*len(V) for _ in range(len(U))]
        D=[[-self.__cost.infini]*len(V) for _ in range(len(U))]
        
        for i in range(1,len(U)):
            D[i][0] = self.__cost.gap_open + (self.__cost.gap_ext * i)
     
        for j in range(1,len(V)):
            I[0][j] = self.__cost.gap_open + (self.__cost.gap_ext * j)
        
        #initialisation des indexs et max score
        max_index=(0,0)
        max_score=-self.__cost.infini
        
        
        for i in range(1,len(U)):
            for j in range(1,len(V)):

                T[i][j] = max(self.__cost(U[i],V[j]) + max(T[i-1][j-1],I[i-1][j-1],D[i-1][j-1]),0)
                
                i1 = self.__cost.gap_open + self.__cost.gap_ext + T[i][j-1]
                i2 = self.__cost.gap_ext + I[i][j-1]
                i3 = self.__cost.gap_open + self.__cost.gap_ext + D[i][j-1]
                I[i][j] = max(i1,i2,i3)
                
                d1 = self.__cost.gap_open + self.__cost.gap_ext + T[i-1][j]
                d2 = self.__cost.gap_open + self.__cost.gap_ext + I[i-1][j]
                d3 = self.__cost.gap_ext + D[i-1][j]
                D[i][j] = max(d1,d2,d3)
                
                if T[i][j]>= max_score :
                    max_score = T[i][j]
                    max_index = (i,j)
         
#        print(pd.DataFrame(np.array(T).transpose(), index = [*V], columns = [*U]))
#        print(pd.DataFrame(np.array(I).transpose(), index = [*V], columns = [*U]))
#        print(pd.DataFrame(np.array(D).transpose(), index = [*V], columns = [*U]))
#        print(max_score,max_index)
        return T,I,D,max_index,max_score
    
    def __align_traceback_gap(self):
        U=" "+ self.__SeqU.seq
        V=" "+ self.__SeqV.seq
        out=[]
        stack=[self.__max_index_T]
        
        while len(stack)>0:
            i,j=stack.pop()
            if i>0 and j>0:
                if self.__matrix[i][j]==self.__cost(U[i],V[j]) + self.__matrix[i-1][j-1]:
                    stack.append((i-1,j-1))
                    if self.__cost(U[i],V[j])==self.__cost.identity:
                        out=[(U[i],':',V[j])]+out
                        self.__aln_stats['identity']+=1
                    else:
                        if self.__cost.intersect_IUPAC(U[i],V[j])>0:
                            out=[(U[i],'.',V[j])]+out
                            self.__aln_stats['similarity']+=1
                        else:
                            out=[(U[i],' ',V[j])]+out
                            
                elif self.__matrix[i][j]==self.__cost(U[i],V[j]) + self.__matrix_I[i-1][j-1]:
                    stack.append((i,j-1))
                    out=[('-',' ',V[j])]+out
                    self.__aln_stats['gaps']+=1
                    
                elif self.__matrix[i][j]==self.__cost(U[i],V[j]) + self.__matrix_D[i-1][j-1]:
                    stack.append((i-1,j))
                    out=[(U[i],' ','-')]+out
                    self.__aln_stats['gaps']+=1
                    
        self.__aln_stats['length']=len(out)
        self.__aln_stats['similarity']+=self.__aln_stats['identity']
        return out,(i,j)
    
    def __align_traceback_rec(self,i,j,out):
        U=" "+ self.__SeqU.seq
        V=" "+ self.__SeqV.seq
        if i>0 and j>0 and self.__matrix[i][j]>0:
            if self.__matrix[i][j]==self.__matrix[i-1][j-1]+self.__cost(U[i],V[j]):
                if self.__cost(U[i],V[j])==self.__cost.identity:
                    outval=[(U[i],':',V[j])]+out
                    self.__aln_stats['identity']+=1
                else:
                    if self.__cost.intersect_IUPAC(U[i],V[j])>0:
                        outval=[(U[i],'.',V[j])]+out
                        self.__aln_stats['similarity']+=1
                    else:
                        outval=[(U[i],' ',V[j])]+out
                return self.__align_traceback_rec(i-1,j-1,outval)

            elif self.__matrix[i][j]==self.__matrix[i-1][j]+self.__cost.indel:
                self.__aln_stats['gaps']+=1
                return self.__align_traceback_rec(i-1,j,[(U[i],' ','-')]+out)
                
            elif self.__matrix[i][j]==self.__matrix[i][j-1]+self.__cost.indel:
                self.__aln_stats['gaps']+=1
                return self.__align_traceback_rec(i,j-1,out=[('-',' ',V[j])]+out)
            
        self.__aln_stats['length']=len(out)
        self.__aln_stats['similarity']+=self.__aln_stats['identity']
        return out,(i,j)

    def __align_traceback(self):
        U=" "+ self.__SeqU.seq
        V=" "+ self.__SeqV.seq
        out=[]
        i,j=self.__max_index
        
        while i>0 and j>0 and self.__matrix[i][j]>0:
            if self.__matrix[i][j]==self.__matrix[i-1][j-1]+self.__cost(U[i],V[j]):
                if self.__cost(U[i],V[j])==self.__cost.identity:
                    out=[(U[i],':',V[j])]+out
                    self.__aln_stats['identity']+=1
                else:
                    if self.__cost.intersect_IUPAC(U[i],V[j])>0:
                        out=[(U[i],'.',V[j])]+out
                        self.__aln_stats['similarity']+=1
                    else:
                        out=[(U[i],' ',V[j])]+out
                i=i-1
                j=j-1
            elif self.__matrix[i][j]==self.__matrix[i-1][j]+self.__cost.indel:
                out=[(U[i],' ','-')]+out
                self.__aln_stats['gaps']+=1
                i=i-1
            elif self.__matrix[i][j]==self.__matrix[i][j-1]+self.__cost.indel:
                out=[('-',' ',V[j])]+out
                self.__aln_stats['gaps']+=1
                j=j-1
                
        self.__aln_stats['length']=len(out)
        self.__aln_stats['similarity']+=self.__aln_stats['identity']
        return out,(i,j)
    
    def __align_traceback_stack(self):
        U=" "+ self.__SeqU.seq
        V=" "+ self.__SeqV.seq
        out=[]
        stack=[self.__max_index]
        
        while len(stack)>0:
            i,j=stack.pop()
            if i>0 and j>0:
                if self.__matrix[i][j]==self.__matrix[i-1][j-1]+self.__cost(U[i],V[j]):
                    stack.append((i-1,j-1))
                    if self.__cost(U[i],V[j])==self.__cost.identity:
                        out=[(U[i],':',V[j])]+out
                        self.__aln_stats['identity']+=1
                    else:
                        if self.__cost.intersect_IUPAC(U[i],V[j])>0:
                            out=[(U[i],'.',V[j])]+out
                            self.__aln_stats['similarity']+=1
                        else:
                            out=[(U[i],' ',V[j])]+out
                elif self.__matrix[i][j]==self.__matrix[i-1][j]+self.__cost.indel:
                    stack.append((i-1,j))
                    out=[(U[i],' ','-')]+out
                    self.__aln_stats['gaps']+=1
                elif self.__matrix[i][j]==self.__matrix[i][j-1]+self.__cost.indel:
                    stack.append((i,j-1))
                    out=[('-',' ',V[j])]+out
                    self.__aln_stats['gaps']+=1
                    
        self.__aln_stats['length']=len(out)
        self.__aln_stats['similarity']+=self.__aln_stats['identity']
        return out,(i,j)
    
    def __print_Aln(self,outAln,min_index,lineSize,outfile=""):
        
        #{'identity':0,'similarity':0,'gaps':0,'length':0}
        identity = self.__aln_stats['identity']
        similarity = self.__aln_stats['similarity']
        gaps = self.__aln_stats['gaps']
        length = self.__aln_stats['length']
            
        out_alignment=""        
        out_alignment+="########################################\n"
        out_alignment+="# Program: {} - v{}\n".format(__title__,__version__)
        out_alignment+="# Rundate: {}\n".format(datetime.now())
        out_alignment+="# Report_file: {report_file}\n"
        out_alignment+="########################################\n"
        out_alignment+="#=======================================\n"
        out_alignment+="# Aligned_sequences: 2\n"
        out_alignment+="# 1: {} - length={}\n".format(self.__SeqU.id,len(self.__SeqU.seq))
        out_alignment+="# 2: {} - length={}\n".format(self.__SeqV.id,len(self.__SeqV.seq))
        out_alignment+="# Gap Affine: {}\n".format(self.__GapAffine)
        out_alignment+="# Costs:\n"
        out_alignment+="#   Indentity: {}\n".format(self.__cost.identity)
        out_alignment+="#   Substitution: {}\n".format(self.__cost.substitution)
        out_alignment+="#   Indel: {}\n".format(self.__cost.indel)
        out_alignment+="#   Gap_penalty: {}\n".format(self.__cost.gap_open)
        out_alignment+="#   Extend_penalty: {}\n".format(self.__cost.gap_ext)
        out_alignment+="#\n"
        out_alignment+="# Length: {}\n".format(self.__aln_stats['length'])
        out_alignment+="# Identity: {}/{} ({:2.1f}%)\n".format(identity,length,float(identity/length)*100)
        out_alignment+="# Similarity: {}/{} ({:2.1f}%)\n".format(similarity,length,float(similarity/length)*100)
        out_alignment+="# Gaps: {}/{} ({:2.1f}%)\n".format(gaps,length,float(gaps/length)*100)
        out_alignment+="# Score: {}\n".format(self.__aln_score)
        out_alignment+="#\n"
        out_alignment+="#=======================================\n\n"
        
        out_alignment+=self.__MultipleLineAln(outAln,min_index,lineSize)
        
        out_alignment+="#---------------------------------------\n"
        out_alignment+="#---------------------------------------\n"

        return out_alignment
    
    def __MultipleLineAln(self,outAln,min_index,lineSize):
        out_alignment=""
        #récupère positions séquences
        #3 pour 3' et 5 pour 5'
        seqU_b3=min_index[0]+1
        seqV_b3=min_index[1]+1
        
        #Copie des séquences pour affichage
        l1Seq=self.__SeqU_Aligned
        l2Info=self.__SeqUV_Matchs
        l3Seq=self.__SeqV_Aligned
        
        maxWidthId=max((len(self.__SeqU.id),len(self.__SeqV.id)))

        while len(l1Seq)>0:
            seqU_b5=seqU_b3+len(l1Seq[:lineSize].replace("-",""))-1
            seqV_b5=seqV_b3+len(l3Seq[:lineSize].replace("-",""))-1
            maxwidth=max([len(str(v)) for v in [seqU_b5,seqU_b3,seqV_b5,seqV_b3]])
            
            out_alignment+="{id:>{widthID}}  {b3:>{width}d}  {seq}  {b5:>{width}d}\n".format(widthID=maxWidthId,width=maxwidth,id=self.__SeqU.id,seq=l1Seq[:lineSize],b3=seqU_b3,b5=seqU_b5)
            out_alignment+="{id:>{widthID}}  {b3:>{width}s}  {seq}  {b5:>{width}s}\n".format(widthID=maxWidthId,width=maxwidth,id=" ",seq=l2Info[:lineSize],b3=" ",b5=" ")
            out_alignment+="{id:>{widthID}}  {b3:>{width}d}  {seq}  {b5:>{width}d}\n".format(widthID=maxWidthId,width=maxwidth,id=self.__SeqV.id,seq=l3Seq[:lineSize],b3=seqV_b3,b5=seqV_b5)
            out_alignment+="\n"
            l1Seq=l1Seq[lineSize:]
            l2Info=l2Info[lineSize:]
            l3Seq=l3Seq[lineSize:]
            seqU_b3=seqU_b5+1
            seqV_b3=seqV_b5+1
            
        return out_alignment
    
    
    #===================
    #Getters Setters
    #===================
    
    @property
    def SeqU_Aligned(self):
        return self.__SeqU_Aligned
    
    @property
    def SeqV_Aligned(self):
        return self.__SeqV_Aligned
    
    @property
    def SeqUV_Matchs(self):
        return self.__SeqUV_Matchs
    
    @property
    def aln_score(self):
        return self.__aln_score
    
    @property
    def SeqU(self):
        return self.__SeqU
    
    @property
    def SeqV(self):
        return self.__SeqV
    
    @property
    def GapAffine(self):
        return self.__GapAffine
    
    @property
    def cost(self):
        return self.__cost
    
    #===================
    #Méthodes magiques
    #=================== 
    
    def __str__(self):
        return self.__MultipleLineAln(self.__outAln,self.__min_index,self.__lineSize)