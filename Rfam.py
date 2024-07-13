#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 09:30:35 2022

@author: vitor
"""

import pandas as pd
import pickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

 
dici2={}
path = os.path.realpath(os.path.dirname(__file__))
rfam1 = pd.read_csv(path+'/family.txt', sep='\t',encoding='cp1252', names=['ID1', 'FamilyName', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27','28', '29', '30', '31', '32', '33'])
rfam1 = rfam1[['ID1', 'FamilyName','17']]

rfam1['17']= rfam1['17'].replace('Gene; rRNA;', 'rRNA')
rfam1['17']= rfam1['17'].replace('Gene; snRNA; snoRNA; scaRNA;', 'snoRNA')
rfam1['17']= rfam1['17'].replace('Gene; snRNA; snoRNA; HACA-box;', 'snoRNA')
rfam1['17']= rfam1['17'].replace('Gene; snRNA; snoRNA; CD-box;', 'snoRNA')
rfam1['17']= rfam1['17'].replace('Gene; snRNA; snoRNA;', 'snoRNA')
rfam1['17']= rfam1['17'].replace('Gene; snRNA; splicing;', 'snRNA')
rfam1['17']= rfam1['17'].replace('Gene; tRNA;', 'tRNA')
rfam1['17']= rfam1['17'].replace('Gene; ribozyme;', 'ribozyme')
rfam1['17']= rfam1['17'].replace('Gene; antisense; ', 'antisense')
rfam1['17']= rfam1['17'].replace('Gene; antisense;', 'antisense')
rfam1['17']= rfam1['17'].replace('Gene; lncRNA;', 'lncRNA')
rfam1['17']= rfam1['17'].replace('Gene; antitoxin;', 'antitoxin')
rfam1['17']= rfam1['17'].replace('Gene; miRNA;', 'miRNA')
rfam1['17']= rfam1['17'].replace('Gene; snRNA;', 'snRNA')
rfam1['17']= rfam1['17'].replace('Gene; sRNA; ', 'sRNA')
rfam1['17']= rfam1['17'].replace('Gene; CRISPR;', 'CRISPR')
rfam1['17']= rfam1['17'].replace('Cis-reg; leader;', 'leader')
rfam1['17']= rfam1['17'].replace('Cis-reg; IRES;', 'IRES')
rfam1['17']= rfam1['17'].replace('Cis-reg; riboswitch;', 'riboswitch')
rfam1['17']= rfam1['17'].replace('Cis-reg; thermoregulator;','thermoregulator')
rfam1['17']= rfam1['17'].replace('Cis-reg; thermoregulator; ','thermoregulator')
rfam1['17']= rfam1['17'].replace('Cis-reg; frameshift_element;', 'frameshift_element')
rfam1['17']= rfam1['17'].replace('Gene;','Gene')
rfam1['17']= rfam1['17'].replace( 'Intron;','Intron')
rfam1['17']= rfam1['17'].replace('Cis-reg; ', 'Cis-reg')
rfam1['17']= rfam1['17'].replace('Cis-reg;', 'Cis-reg')
rfam1['17']= rfam1['17'].replace('Gene; sRNA;', 'sRNA')    
   
rfam2 = rfam1[['ID1', 'FamilyName', '17']]
rfam2.to_csv(path+'/Rfam_dict.csv', index=False)
rfam3 = rfam2[['ID1', '17']]
    
dici2 = {}
dici2 = dict(rfam3.values)
    
with open(path+"/Rfam_dicionario.pkl", "wb") as tf:
    pickle.dump(dici2,tf)

