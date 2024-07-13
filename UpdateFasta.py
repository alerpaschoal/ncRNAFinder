#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 13:27:09 2022

@author: vitor
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import pandas as pd
import pickle

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--fasta", required=True, help= "Original FASTA file")
args = vars(ap.parse_args())

frags=[]
ids=[]
n=1
for record in SeqIO.parse(args['fasta'], 'fasta'):
    ids.append(record.id)
    frag = record.seq
    gravar = SeqRecord(frag, 'Seq_'+str(n), "", "")
    frags.append(gravar)
    n=n+1
SeqIO.write(frags, 'New_'+args['fasta'], 'fasta')

seqDmel=[]
for i in range(len(ids)):
    seqDmel.append('Seq_'+str(i+1))
dfDmel = pd.DataFrame(columns=['Seq','Original'])
dfDmel.Original=ids
dfDmel.Seq = seqDmel

dici=dict(dfDmel.values)

with open("FASTA_Temp.pkl", "wb") as tf:
    pickle.dump(dici,tf)

