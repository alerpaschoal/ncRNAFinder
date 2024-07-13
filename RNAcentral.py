#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 23:27:00 2022

@author: vitor
"""

import pandas as pd
import pickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import gzip

dici={}
path = os.path.realpath(os.path.dirname(__file__))
rcen = pd.read_csv(path+'/id_mapping.tsv', sep='\t', names=['ID', 'ENA', 'GU', 'FimID', 'Type', 'nada'])
rcen = rcen[['ID', 'FimID','Type']]
    
rcen['FimID'] = rcen.FimID.astype('str')
rcen["Id"] = rcen['ID'].map(str) + "_" + rcen["FimID"]
    
rcen = rcen[['Id', 'Type']]
print(rcen.head())

dici=dict(rcen.values) 
    
with open(path+"/RNAcentralV23_ncRNAs_specific.pkl", "wb") as tf:
    pickle.dump(dici,tf)
dici={}
    
   
dict2={}


for seq_record in SeqIO.parse(path+'/rnacentral_active.fasta', "fasta"):
    frase = seq_record.description.split()
    dict2[frase[0]]=frase[1]

with open(path+"/RNAcentralV23_active.pkl", "wb") as tf2:
    pickle.dump(dict2,tf2) 

dict3={}

for seq_record in SeqIO.parse(path+'/rnacentral_active.fasta', "fasta"):
    frase = seq_record.description.split()
    compr = len(seq_record.seq)
    dict3[frase[0]]=compr
    print(dict3)

with open(path+"/RNAcentralV24_active_len.pkl", "wb") as tf3:
    pickle.dump(dict3,tf3) 

print('fim')
