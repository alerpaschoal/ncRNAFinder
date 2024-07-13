#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27/06/2024

@author: vitor

Otimizando execução das funções 
"""

import pandas as pd
import datetime as dt
import argparse
import os
import sys
import pickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
pd.set_option('mode.chained_assignment', None)
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn2, venn2_unweighted
import csv
from joblib import Parallel, delayed 
import time


ap = argparse.ArgumentParser()
ap.add_argument("-b", "--blast", required=True, help= "Result of blast")
ap.add_argument("-i", "--infernal", required=True, help= "Result of infernal")
ap.add_argument("-f", "--fasta", required=False, help= "Original FASTA file")
ap.add_argument("-c", "--cover", type=int, required=False, choices=range(0,101), metavar="[0-100]", default= "95", help= "Coverage value")
ap.add_argument("-r", "--BestHit", required=False, choices=['1', '0'], default= "1", help= "Select the Best result for diferent strands: 1 = Yes ou 0 = No")
ap.add_argument("-t", "--threads", required=False, default= "1", help= "Number of Threads")


args = vars(ap.parse_args())


##############################################################INICIO DAS FUNÇÕES####################################################################
def overlaps_blast(arq):
    remover = set()
    #vet = arq['targetName'].unique()
    
    # Iterar sobre cada targetName único
    #for x in vet:
    teste = arq
    seqFrom = teste['seqFrom'].values
    seqTo = teste['seqTo'].values
    pident = teste['pident'].values
    qcover = teste['qcover'].values
    Evalue = teste['Evalue'].values
    indices = teste.index.values
    
    # Usar arrays NumPy para comparação vetorizada
    for i in range(len(indices)):
        if indices[i] in remover:
            continue
        overlap_mask = (
            (seqFrom[i] >= seqFrom) & (seqFrom[i] <= seqTo) &
            (seqTo[i] >= seqFrom) & (seqTo[i] <= seqTo)
        ) | (
            (seqFrom[i] >= seqFrom) & (seqFrom[i] <= seqTo) &
            (seqTo[i] >= seqFrom) & (seqTo[i] >= seqTo)
        )
        
        candidates = np.where(overlap_mask)[0]
        for j in candidates:
            if i != j:
                if pident[i] == pident[j]:
                    if qcover[i] == qcover[j]:
                        if Evalue[i] > Evalue[j]:
                            remover.add(indices[i])
                        else:
                            remover.add(indices[j])
                    elif qcover[i] > qcover[j]:
                        remover.add(indices[j])
                    else:
                        remover.add(indices[i])
                elif pident[i] > pident[j]:
                    remover.add(indices[j])
                else:
                    remover.add(indices[i])

    return arq.drop(remover)



def overlaps(arq):
    remover=set()
    teste = arq

    seqFrom = teste['seqFrom'].values
    seqTo = teste['seqTo'].values
    score = teste['score'].values
    indices = teste.index.values

    for i in range(len(indices)):
        if indices[i] in remover:
            continue
        overlap_mask = (
            (seqFrom[i] >= seqFrom) & (seqFrom[i] <= seqTo) &
            (seqTo[i] >= seqFrom) & (seqTo[i] <= seqTo)
        ) | (
            (seqFrom[i] >= seqFrom) & (seqFrom[i] <= seqTo) &
            (seqTo[i] >= seqFrom) & (seqTo[i] >= seqTo)
        )

        candidates = np.where(overlap_mask)[0]
        for j in candidates:
                if i != j:
                    if score[i] > score[j]:
                        remover.add(indices[j])
                    else:
                        remover.add(indices[i])

    return arq.drop(remover)


def overlaps_blast_infernal(arqBlast, arqInfernal):
    operacao = 'U'
    remover = []

    def calcular_intervalos(seqFrom1, seqTo1, seqFrom2, seqTo2):
        maiorFrom = max(seqFrom1, seqFrom2)
        menorFrom = min(seqFrom1, seqFrom2)
        maiorTo = max(seqTo1, seqTo2)
        menorTo = min(seqTo1, seqTo2)
        return maiorFrom, menorFrom, maiorTo, menorTo

    def checar_sobreposicao(seqFrom1, seqTo1, seqFrom2, seqTo2):
        return seqFrom1 <= seqTo2 and seqTo1 >= seqFrom2

    for i in arqBlast.index:
        for j in arqInfernal.index:
            if arqBlast.loc[i, 'targetName'] == arqInfernal.loc[j, 'targetName']:
                if checar_sobreposicao(arqBlast.loc[i, 'seqFrom'], arqBlast.loc[i, 'seqTo'], arqInfernal.loc[j, 'seqFrom'], arqInfernal.loc[j, 'seqTo']):
                    maiorFrom, menorFrom, maiorTo, menorTo = calcular_intervalos(arqBlast.loc[i, 'seqFrom'], arqBlast.loc[i, 'seqTo'], arqInfernal.loc[j, 'seqFrom'], arqInfernal.loc[j, 'seqTo'])

                    if arqBlast.loc[i, 'queryName'] == arqInfernal.loc[j, 'queryName'] and arqBlast.loc[i, 'source'] != arqInfernal.loc[j, 'source']:
                        remover.append(j)
                        arqBlast.loc[i, 'source'] = 'Both'
                    elif arqBlast.loc[i, 'queryName'] != arqInfernal.loc[j, 'queryName'] and arqBlast.loc[i, 'source'] != arqInfernal.loc[j, 'source']:
                        remover.append(j)
                        arqBlast.loc[i, 'source'] = 'Both'
                        arqBlast.loc[i, 'ncRNA'] = 'ncRNA'

                    if operacao == "I":
                        arqBlast.loc[i, 'seqFrom'] = maiorFrom
                        arqBlast.loc[i, 'seqTo'] = menorTo
                    else:
                        arqBlast.loc[i, 'seqFrom'] = menorFrom
                        arqBlast.loc[i, 'seqTo'] = maiorTo    
    return arqBlast, arqInfernal.drop(remover)



def overlaps_BestHit(arq):
    remover=set()
    
    teste = arq
    seqFrom = teste['seqFrom'].values
    seqTo = teste['seqTo'].values
    Score = teste['score'].values
    Evalue = teste['Evalue'].values
    indices = teste.index.values

    for i in range(len(indices)):
        if indices[i] in remover:
            continue
        overlap_mask = (
            (seqFrom[i] >= seqFrom) & (seqFrom[i] <= seqTo) &
            (seqTo[i] >= seqFrom) & (seqTo[i] <= seqTo)
        ) | (
            (seqFrom[i] >= seqFrom) & (seqFrom[i] <= seqTo) &
            (seqTo[i] >= seqFrom) & (seqTo[i] >= seqTo)
        )
        
        candidates = np.where(overlap_mask)[0]
        for j in candidates:
            if i != j:
                if Evalue[i] > Evalue[j]:
                    remover.add(indices[i])
                elif Evalue[i] < Evalue[j]:
                    remover.add(indices[j])
                else:
                    if Score[i] >= Score[j]:
                        remover.add(indices[j])
                    else:
                        remover.add(indices[i])

    return arq.drop(remover)

###Salvar tabela com resultados originais
def salvar_original(resultFinal, arqBlast1, arqInfernal):
    arqBlast1 = arqBlast1.replace({'minus': '-', 'plus':'+'})
    for index6, row6 in arqBlast1.iterrows():
        if arqBlast1.loc[index6,'seqFrom'] > arqBlast1.loc[index6,'seqTo']:
            aux6 = arqBlast1.loc[index6,'seqFrom']
            arqBlast1.loc[index6,'seqFrom'] = arqBlast1.loc[index6,'seqTo'] 
            arqBlast1.loc[index6,'seqTo'] = aux6
    
    if os.path.exists(nameI+"_ID_annotation.txt"):
        os.remove(nameI+"_ID_annotation.txt")
    
    for j in resultFinal.index:
        if resultFinal.loc[j,'source'] == 'Blast':
            sTarget = resultFinal['targetName'][j]
            sFrom = resultFinal['seqFrom'][j]
            sTo = resultFinal['seqTo'][j]
            blast = arqBlast1.query('seqFrom>=@sFrom & seqTo<=@sTo & targetName==@sTarget')
            blast2 = blast.rename(columns={'IDdb':'qseqid', 'targetName':'sseqid', 'pident':'pident', 'len':'lenght', 'mdl': 'mismatch', 'mdlFrom':'gapopen', 'mdlTo':'qstart', 'trunc':'qend', 'score':'bitscore', 'inc':'qlen', 'bias':'slen'})
            
            salvar = open(nameI+'_ID_annotation.txt', 'a+')
            salvar.write(resultFinal.loc[j,'ID'] +'\n'+'Blast: '+'\n')
            salvar.close()
            salvar = open(nameI+'_ID_annotation.txt', 'a+')
            blast2.to_csv(nameI+'_ID_annotation.txt', sep='\t',mode='a+', index=False, header=False)
            salvar.write('\n')
            salvar.close()
        elif resultFinal.loc[j,'source'] == 'Infernal':
            sTarget = resultFinal['targetName'][j]
            sFrom = resultFinal['seqFrom'][j]
            sTo = resultFinal['seqTo'][j]
            Infernal = arqInfernal.query('seqFrom>=@sFrom & seqTo<=@sTo & targetName==@sTarget')
            Infernal = Infernal[['targetName', 'seqFrom', 'seqTo', 'pident', 'qcover', 'IDdb', 'strand', 'score', 'Evalue', 'ncRNA']]
            
            salvar = open(nameI+'_ID_annotation.txt', 'a+')
            salvar.write(resultFinal.loc[j,'ID'] +'\n'+'Infernal: '+'\n')
            salvar.close()
            salvar = open(nameI+'_ID_annotation.txt', 'a+')
            Infernal.to_csv(nameI+'_ID_annotation.txt', sep='\t',mode='a+', index=False, header=False)
            salvar.write('\n')
            salvar.close()
        else:
            sTarget = resultFinal['targetName'][j]
            sFrom = resultFinal['seqFrom'][j]
            sTo = resultFinal['seqTo'][j]
            blast = arqBlast1.query('seqFrom>=@sFrom & seqTo<=@sTo & targetName==@sTarget')
            
            Infernal = arqInfernal.query('seqFrom>=@sFrom & seqTo<=@sTo & targetName==@sTarget')
            blast2 = blast.rename(columns={'IDdb':'qseqid', 'targetName':'sseqid', 'pident':'pident', 'queryName':'lenght', 'mdl': 'mismatch', 'mdlFrom':'gapopen', 'mdlTo':'qstart', 'trunc':'qend', 'score':'bitscore', 'inc':'qlen', 'bias':'slen'})
            for index8, row8 in Infernal.iterrows():
                if Infernal.loc[index8,'strand'] == '-' or Infernal.loc[index8,'strand'] == 'minus':
                    if Infernal.loc[index8,'seqFrom'] < Infernal.loc[index8,'seqTo']:
                        aux8 = Infernal.loc[index8,'seqFrom']
                        Infernal.loc[index8,'seqFrom'] = Infernal.loc[index8,'seqTo'] 
                        Infernal.loc[index8,'seqTo'] = aux8        
            Infernal = Infernal.drop_duplicates(subset=['targetName', 'seqTo', 'seqFrom'])
                                   
            salvar = open(nameI+'_ID_annotation.txt', 'a+')
            salvar.write(resultFinal.loc[j,'ID'] +'\n'+'Blast: '+'\n')
            salvar.close()
            salvar = open(nameI+'_ID_annotation.txt', 'a+')
            blast2.to_csv(nameI+'_ID_annotation.txt', sep='\t',mode='a+', index=False, header=False)
            salvar.write('Infernal: '+'\n')
            salvar.close()
            salvar = open(nameI+'_ID_annotation.txt', 'a+')
            salvar.write('\n')
            Infernal.to_csv(nameI+'_ID_annotation.txt', sep='\t',mode='a+', index=False, header=False)
            salvar.close()


###Salva em GFF        
def salvar_gff(resultFinalFinal, nameI):      
    if os.path.exists(nameI+'_annotation.gff'):
        os.remove(nameI+'_annotation.gff')
    
    gff = pd.DataFrame(columns=['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
    gff['seq_id'] = resultFinalFinal['targetName']
    gff['source'] = resultFinalFinal['source']
    gff['type'] = resultFinalFinal['queryName']
    gff['start'] = resultFinalFinal['seqFrom']
    gff['end'] = resultFinalFinal['seqTo']
    gff['score'] = resultFinalFinal['score']
    gff['strand'] = resultFinalFinal['strand']
    gff['phase'] = '.'
    gff['attributes'] = resultFinalFinal['ID'] #resultado do blast/infernal
    gff['attributes']=gff['attributes'].str.replace('id', 'ID=id')
    gff.to_csv(nameI+'_annotation.gff', sep='\t',mode='a', index=False, header=False)
    
    file = open(nameI+'_annotation.gff', 'r')
    lines = file.readlines()
    file.close()
    lines.insert(0, "##gff-version 3"+'\n')
    
    file = open(nameI+"_annotation.gff", 'w')
    file.writelines(lines)
    file.close()
    
    return gff


def salvar_fasta(gff, nameI):
    ids = list(gff['seq_id'])
    start = list(gff['start'])
    end = list(gff['end'])
    attr = list(gff['attributes'].str.replace('ID=id',''))
    attrOri = list(gff['attributes'])
    tipo = list(gff['type'])
    frags=[]
    for record in SeqIO.parse(args['fasta'], "fasta"):
        for i in range(len(ids)):
            if record.id == ids[i]:
                frag = record.seq[start[i]-1:end[i]]
                gravar = SeqRecord(frag, "fragment_" + str(attr[i]), ids[i], tipo[i]+';'+attrOri[i])
                frags.append(gravar)
    SeqIO.write(frags, nameI+"_sequences_ncRNA.fa", "fasta")

def qtdTotal_ncRNA(resultFinal, nameI):
    QtdTotalRNA = resultFinal['queryName'].value_counts(dropna=False, sort=True)
    QtdTotalRNA = QtdTotalRNA.to_frame().reset_index()
    QtdTotalRNA = QtdTotalRNA.set_axis(['ncRNA', 'Qtd'], axis=1)
    QtdTotalRNA.to_csv(nameI+'_Table_quantity_ncRNAs.csv')
    
def graficoTotal_barra(resultFinal, nameI):
    plt.clf()
    QtdTotalRNA = resultFinal['queryName'].value_counts(dropna=False, sort=True)
    QtdTotalRNA = QtdTotalRNA.to_frame().reset_index()
    QtdTotalRNA = QtdTotalRNA.set_axis(['ncRNA', 'Qtd'], axis=1)
        
    fig, ax = plt.subplots()
    my_cmap = plt.get_cmap("Set1")
    graf = ax.bar(QtdTotalRNA['ncRNA'].to_list(), QtdTotalRNA['Qtd'].to_list(), color=my_cmap.colors)
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Quantidade')
    ax.set_title('Quantidade x ncRNA - '+ nameI)
    my_cmap = plt.get_cmap("viridis")
    ax.set_xticks(np.arange(len(QtdTotalRNA['ncRNA'])), QtdTotalRNA['ncRNA'], rotation=45)
    ax.bar_label(graf, padding=1)
    fig.tight_layout()
    plt.savefig(nameI+'_BarPlot.png', format='png', dpi=1200)
    

def graficoTotal_stacked(resultFinal, nameI):
    plt.clf()
    vet = resultFinal['targetName'].unique()
    num=0
    data = pd.DataFrame()
    for i in vet:
        teste1 = resultFinal[resultFinal['targetName']==i]
        if num==0:
            aux = teste1['queryName'].value_counts(dropna=False, sort=True)
            aux = aux.to_frame(name = i)
            data = aux
            num=1
        else:
            aux1 = teste1['queryName'].value_counts(dropna=False, sort=True)
            aux1 = aux1.to_frame(name = i)
            data = pd.merge(data, aux1, how='outer', right_index = True, left_index = True)
            
    data.fillna(value = 0, inplace = True)
    
    
    data.T.plot(kind="bar", colormap="tab20", stacked=True,figsize=(10,8), title = nameI).get_legend().set_visible(True) #não está plotando com matplotlib
    plt.savefig(nameI+'_StackedPlot.png', format='png', dpi=1200)

def diagrama_venn(resultFinal, nameI):
    plt.clf()
    QtdFerramenta = resultFinal['source'].value_counts(dropna=False, sort=True)
    QtdFerramenta = QtdFerramenta.to_frame().T
    
    try:
        set3 = int(QtdFerramenta['Both'])
    except:
        set3=0
    try:
        set2 = int(QtdFerramenta['Infernal'])
    except:
        set2=0
    try:
        set1 = int(QtdFerramenta['Blast'])
    except:
        set1=0
    
    plt.title(nameI)
    out = venn2_unweighted(subsets = (set1, set2, set3), set_labels = ('Blast', 'Infernal'), set_colors=('purple', 'skyblue'))
    for text in out.set_labels:
       text.set_fontsize(25)
    for text in out.subset_labels:
       text.set_fontsize(18)
    plt.savefig(nameI+'_DiagramaVeen.png', format='png', dpi=1200)

def resultados_tRNA(resultFinal, nameI):
    tRNA = resultFinal[resultFinal['queryName']=="tRNA"]
    tRNA = tRNA.reset_index(drop=True)
    if tRNA.empty == False:
        tRNA.to_csv(nameI+'_tRNAs.csv', index=False)
        gfftRNA = salvar_gff(tRNA, nameI+'_tRNA')
        QtdTotaltRNA = tRNA['ncRNA'].value_counts(dropna=False, sort=True)
        QtdTotaltRNA = QtdTotaltRNA.to_frame().reset_index()
        QtdTotaltRNA = QtdTotaltRNA.set_axis(['ncRNA', 'Qtd'], axis=1)
        
        tRNA2 = tRNA[['ncRNA', 'IDdb']]
        tRNA2 = tRNA2.drop_duplicates()
        tRNA2 = tRNA2.set_index('ncRNA')
        tRNAdict = tRNA2.to_dict()
        QtdTotaltRNA['IDdb'] = QtdTotaltRNA['ncRNA'].map(tRNAdict['IDdb'])
        
        QtdTotaltRNA.to_csv(nameI+'_Table_quantity_tRNAs.csv')
        if args['fasta'] != None:
            salvar_fasta(gfftRNA, nameI+'_tRNA')
        
        
    
def resultados_miRNA(resultFinal, nameI):
    miRNA = resultFinal[(resultFinal['queryName']=="miRNA") | (resultFinal['queryName']=="pre-miRNA")]
    miRNA = miRNA.reset_index(drop=True)
    if miRNA.empty == False:
        miRNA.to_csv(nameI+'_miRNAs.csv', index=False)
        gffmiRNA = salvar_gff(miRNA, nameI+'_miRNA')
        QtdTotalmiRNA = miRNA['ncRNA'].value_counts(dropna=False, sort=True)
        QtdTotalmiRNA = QtdTotalmiRNA.to_frame().reset_index()
        QtdTotalmiRNA = QtdTotalmiRNA.set_axis(['ncRNA', 'Qtd'], axis=1)
        
        miRNA2 = miRNA[['ncRNA', 'IDdb']]
        miRNA2 = miRNA2.drop_duplicates()
        miRNA2 = miRNA2.set_index('ncRNA')
        miRNAdict = miRNA2.to_dict()
        QtdTotalmiRNA['IDdb'] = QtdTotalmiRNA['ncRNA'].map(miRNAdict['IDdb'])
        
        QtdTotalmiRNA.to_csv(nameI+'_Table_quantity_miRNAs.csv')
        if args['fasta'] != None:    
            salvar_fasta(gffmiRNA, nameI+'_miRNA')
 
        
 
def pipeBlast(arqBlast1):
    
    #arqBlast1['qcover'] = round(100 *((abs(arqBlast1['seqTo'] - arqBlast1['seqFrom']) +1) /arqBlast1['inc']),3)
    arqBlast1['qcover'] = round(100 *((abs(arqBlast1['trunc'] - arqBlast1['mdlTo'] +1)) /arqBlast1['inc']),3)
    cover = int(args['cover'])
    arqBlast1= arqBlast1.query('qcover>= @cover')

    infile = open("FASTA_Temp.pkl",'rb')
    FastaOriginal = pickle.load(infile)
    infile.close()
        
    arqBlast1['targetName'] = arqBlast1['Original'].map(FastaOriginal)

    arqBlast1 = arqBlast1.sort_values(['targetName', 'pident', 'qcover', 'Evalue', 'seqFrom', 'seqTo', 'score'], ascending=[True, False, False, True, True, False, False]) #ordena as sequencias
    arqBlast1 = arqBlast1.drop_duplicates(subset=['targetName', 'seqFrom', 'seqTo', 'strand']) #remove quando tem inicio e fim igual, mantendo a com menor e-value
    arqBlast1 = arqBlast1.drop_duplicates(subset=['targetName', 'seqFrom', 'strand']) #remove quando tem inicio igual, mantendo a com menor e-value
    arqBlast1 = arqBlast1.drop_duplicates(subset=['targetName', 'seqTo', 'strand']) #remove quando tem fim igual, mantendo a com menor e-value
    arqBlast1 = arqBlast1.reset_index(drop=True) #reenumera o index do dataframe

    #arruma as coordenadas de inicio e fim para o menor valor ficar no inicio e o maior no fim
    for index, row in arqBlast1.iterrows():
        if arqBlast1.loc[index,'seqFrom'] > arqBlast1.loc[index,'seqTo']:
            aux = arqBlast1.loc[index,'seqFrom']
            arqBlast1.loc[index,'seqFrom'] = arqBlast1.loc[index,'seqTo'] 
            arqBlast1.loc[index,'seqTo'] = aux
            
   
    #abre o dicionario com os id e seus respectivos ncRNAs
    try:
        if os.path.exists(path+"/RNAcentralV23_ncRNAs_specific.pkl"):
            infile = open(path+"/RNAcentralV23_ncRNAs_specific.pkl",'rb')
            dici1 = pickle.load(infile)
            infile.close()
    except: 
        print('Falta RNAcentral_ncRNAs_specific.pkl')
    
    try:
        if os.path.exists(path+"/RNAcentralV23_active.pkl"):
            infile = open(path+"/RNAcentralV23_active.pkl",'rb')
            dici2 = pickle.load(infile)
            infile.close()
    except: 
        print('Falta RNAcentral_active.pkl')
        
    arqBlast1['ncRNA']=arqBlast1['IDdb'].map(dici1) #adiciona uma coluna com o ncRNA da anotação
    arqBlast1['queryName']=arqBlast1['IDdb'].map(dici2)
    arqBlast = arqBlast1[['targetName', 'seqFrom', 'seqTo', 'pident', 'qcover', 'IDdb', 'strand', 'score', 'Evalue', 'ncRNA', 'queryName']] #copia o dataframe com colunas especificas
    arqBlast = arqBlast.replace({'minus': '-', 'plus':'+'}) #muda para + e - os sentidos das fitas
    
    dici1 ={} #fecha o dicionario
    dici2 ={}
    
    #Separa o dataframe por fita positiva e negativa. E na fita negativa a coluna seqTo é trocada com a seqFrom, pq o bedtools precisa que a primeira coordenada seja menor que a segunda
    negativa = arqBlast[arqBlast['strand']=='-']
    positiva = arqBlast[arqBlast['strand']=='+']
    
    #Salva como resultado do Blast
    negativa['source']= 'Blast'
    positiva['source']= 'Blast'
       
    negativa=negativa.reset_index(drop=True) #arruma o index
    positiva=positiva.reset_index(drop=True) #arruma o index
    
    #resultPositiva = overlaps_blast(positiva) #chama a função para identificar e tratar os overlaps
    #resultNegativa = overlaps_blast(negativa) #chama a função para identificar e tratar os overlaps
    vet1 = positiva['targetName'].unique()
    vet2 = negativa['targetName'].unique()
    resultadoPositiva = Parallel(n_jobs =int(args['threads']))(delayed(overlaps_blast)(positiva[positiva['targetName']==x]) for x in vet1)
    resultadoNegativa = Parallel(n_jobs = int(args['threads']))(delayed(overlaps_blast)(negativa[negativa['targetName']==x]) for x in vet2)

    resultPositiva = pd.concat(resultadoPositiva)
    resultNegativa = pd.concat(resultadoNegativa)

    resultBlast2 = pd.merge(resultPositiva, resultNegativa, how ='outer') #junto os resultados em um dataframe só
    resultBlast2 = resultBlast2[['targetName', 'seqFrom', 'seqTo', 'pident', 'qcover', 'IDdb', 'strand', 'score', 'Evalue', 'ncRNA', 'queryName', 'source']] #organiza o datafraame
    resultBlast2 = resultBlast2.sort_values(['targetName', 'seqFrom', 'seqTo', 'Evalue', 'score'], ascending=[True,True, False, True, False]) #ordena o dataframe
    resultBlast2['ncRNA'] = resultBlast2.ncRNA.astype('object')
    resultBlast2['queryName'] = resultBlast2.queryName.astype('object')
    resultBlast2 = resultBlast2.reset_index(drop=True)  #arruma o index
    
    return resultBlast2, arqBlast1

    
def pipeInf(arq1):
    #arq1 = arq1[:-10] #remove as 10 ultimas linhas
    arq1['qcover'] = '-'
    arq1['seqFrom'] = arq1.seqFrom.astype('int64')
    arq1['seqTo'] = arq1.seqTo.astype('int64')
    arq1['score'] = arq1.score.astype('float64')
    arq1['Evalue'] = arq1.Evalue.astype('float64')
    
    arq1['targetName'] = arq1['Original'] #.map(FastaOriginal)
    
    arq1 = arq1.sort_values(['targetName', 'score', 'seqFrom', 'seqTo', 'Evalue'], ascending=[True,False, True, False, True])
    arq1 = arq1.reset_index(drop=True)
    
    for index2, row2 in arq1.iterrows():
        if arq1.loc[index2,'seqFrom'] > arq1.loc[index2,'seqTo']:
            aux2 = arq1.loc[index2,'seqFrom']
            arq1.loc[index2,'seqFrom'] = arq1.loc[index2,'seqTo'] 
            arq1.loc[index2,'seqTo'] = aux2
        
    arq1 = arq1.drop_duplicates(subset=['targetName', 'seqFrom', 'seqTo'])#remove quando tem inicio e fim igual, mantendo a com menor e-value
    arq1 = arq1.drop_duplicates(subset=['targetName', 'seqFrom'])#remove quando tem inicio igual, mantendo a com menor e-value
    arq1 = arq1.drop_duplicates(subset=['targetName', 'seqTo'])#remove quando tem fim igual, mantendo a com menor e-value
    arqInfernal1 = arq1[['targetName', 'pident', 'queryName', 'IDdb', 'qcover', 'mdlFrom', 'mdlTo', 'seqFrom', 'seqTo', 'strand', 'trunc', 'passs', 'bias','inc', 'score', 'Evalue', 'description']]
    arq = arq1[['targetName', 'seqFrom', 'seqTo', 'pident',  'qcover', 'IDdb', 'strand', 'score', 'Evalue', 'queryName']]
    
    
    arqInfernal = arq.sort_values(['targetName', 'seqFrom', 'seqTo', 'score', 'Evalue'], ascending=[True,True, False, False, True])
    
    try:
        if os.path.exists(path+"/Rfam_dicionario.pkl"):
            infile = open(path+"/Rfam_dicionario.pkl",'rb')
            dici3 = pickle.load(infile)
            infile.close()
    except: 
        print('Falta Rfam_dicionario.pkl')
    
    arqInfernal['ncRNA'] = arqInfernal['queryName']
    arqInfernal['queryName']=arqInfernal['IDdb'].map(dici3)
    arqInfernal['source']= 'Infernal'
    arqInfernal = arqInfernal.sort_values(['targetName', 'score', 'seqFrom', 'seqTo', 'Evalue'], ascending=[True,False, True, False, True])
    arqInfernal = arqInfernal.reset_index(drop=True)
    
    dici3={}
    
    #resultInfernal = overlaps(arqInfernal)
    vet3 = arqInfernal['targetName'].unique()
    resultadoArqInfernal = Parallel(n_jobs = int(args['threads']))(delayed(overlaps)(arqInfernal[arqInfernal['targetName']==x]) for x in vet3)
    resultInfernal = pd.concat(resultadoArqInfernal)

    resultInfernal2 = resultInfernal[['targetName', 'seqFrom', 'seqTo', 'pident',  'qcover', 'IDdb', 'strand', 'score', 'Evalue', 'queryName','ncRNA', 'source']]
    resultInfernal2 = resultInfernal2.sort_values(['targetName', 'seqFrom', 'seqTo', 'Evalue', 'score'], ascending=[True,True, False, False, True])
    #resultInfernal2['ncRNA'] = '-'

    return resultInfernal2, arqInfernal  
        
 

####################################################################FIM DAS FUNÇÕES####################################################################


##Inicia a identificação de overlaps no Blast
path = os.path.realpath(os.path.dirname(__file__))

blast1= args['blast'] #salva o path do arquivo
nameB = blast1.rstrip('.tsv') #salva apenas o inicio do nome do arquivo
#arqBlast1 = pd.read_csv(blast1, sep='\s+', names=['IDdb', 'Original',  'pident', 'qcover', 'len', 'mdl', 'mdlFrom', 'mdlTo', 'trunc','seqFrom', 'seqTo', 'strand', 'Evalue', 'score','inc', 'bias']) #abre o resultado do Blast em um dataframe
chunksize = 100000

# Processar o arquivo em chunks
chunk_list = []  # Lista para armazenar cada chunk processado
for chunk in pd.read_csv(blast1, sep='\s+', names=['IDdb', 'Original',  'pident', 'qcover', 'len', 'mdl', 'mdlFrom', 'mdlTo', 'trunc','seqFrom', 'seqTo', 'strand', 'Evalue', 'score','inc', 'bias'], chunksize=chunksize):
    chunk['qcover'] = round(100 *((abs(chunk['trunc'] - chunk['mdlTo'] +1)) /chunk['inc']),3)
    cover = int(args['cover'])
    chunk= chunk.query('qcover>= @cover')
    
    chunk_list.append(chunk)

# Concatenar todos os chunks em um único DataFrame, se necessário
arqBlast1 = pd.concat(chunk_list, ignore_index=True)


##inicia a analise do Infernal
infernal1= args['infernal'] #path do arquivo do infernal
nameI = infernal1[:-4] #salva o inicio do nome do arquivo
arq1 = pd.read_csv(infernal1, sep=";", quoting=csv.QUOTE_NONE, names=['queryName', 'IDdb', 'Original', 'qcover', 'pident','1','mdlFrom', 'mdlTo', 'seqFrom', 'seqTo', 'strand', 'trunc', 'passs', 'bias','inc', 'score', 'Evalue', 'description', '2', '3', '4', '5', '6', '7', '8', '9'])


if (len(arqBlast1) != 0) & (len(arq1) != 0):
    resultBlast2, arqBlast1 = pipeBlast(arqBlast1)
    resultInfernal2, arqInfernal = pipeInf(arq1)

    BlastFinal, InfernalFinal = overlaps_blast_infernal(resultBlast2, resultInfernal2)
        
    #BlastFinal['seqFrom'] = BlastFinal.seqFrom.astype('int64')
    #BlastFinal['seqTo'] = BlastFinal.seqTo.astype('int64')
    BlastFinal['pident'] = BlastFinal.pident.astype('object')
    BlastFinal['qcover'] = BlastFinal.qcover.astype('object')
        
        
    juntoFinal = pd.merge(BlastFinal, InfernalFinal, how ='outer')

elif (len(arqBlast1) == 0) & (len(arq1) != 0):
    juntoFinal, arqInfernal = pipeInf(arq1)
    arqBlast1=[]
    
elif (len(arqBlast1) != 0 )& (len(arq1) == 0):
    juntoFinal, arqBlast1 = pipeBlast(arqBlast1)
    arqInfernal=[]
    
else:
    print('No result! Exiting...')
    sys.exit()

juntoFinal = juntoFinal[['targetName', 'seqFrom', 'seqTo', 'pident', 'qcover', 'IDdb', 'strand', 'score', 'Evalue','ncRNA','queryName', 'source']]
juntoFinal['seqFrom'] = juntoFinal.seqFrom.astype('int')
juntoFinal['seqTo'] = juntoFinal.seqTo.astype('int')


juntoFinal = juntoFinal.sort_values(['targetName', 'pident', 'qcover', 'Evalue', 'seqFrom', 'seqTo'], ascending=[True,False, False, True, True, False])
juntoFinal = juntoFinal.drop_duplicates(subset=['targetName', 'seqFrom', 'seqTo', 'strand']) #remove quando tem inicio e fim igual, mantendo a com menor e-value
juntoFinal = juntoFinal.drop_duplicates(subset=['targetName', 'seqFrom', 'strand']) #remove quando tem inicio igual, mantendo a com menor e-value
juntoFinal = juntoFinal.drop_duplicates(subset=['targetName', 'seqTo', 'strand']) #remove quando tem fim igual, mantendo a com menor e-value

if args['BestHit'] == '1':
    melhorResultado = juntoFinal.sort_values(['targetName', 'pident', 'qcover', 'Evalue', 'score'], ascending=[True,False, False, True, False])
    melhorResultado = melhorResultado.drop_duplicates(subset=['targetName', 'seqFrom', 'seqTo']) #remove quando tem inicio e fim igual, mantendo a com menor e-value
    melhorResultado = melhorResultado.drop_duplicates(subset=['targetName', 'seqFrom']) #remove quando tem inicio igual, mantendo a com menor e-value
    melhorResultado = melhorResultado.drop_duplicates(subset=['targetName', 'seqTo']) #remove quando tem fim igual, mantendo a com menor e-value
    #melhorResultado = overlaps_BestHit(melhorResultado) #remove quando tiver fitas diferentes, mantendo a com menor e-value ou maior score
    vet4 = melhorResultado['targetName'].unique()
    melhorResultado2 = Parallel(n_jobs = int(args['threads']))(delayed(overlaps_BestHit)(melhorResultado[melhorResultado['targetName']==x]) for x in vet4)
    melhorResultado = pd.concat(melhorResultado2)
    melhorResultado = melhorResultado.sort_values(['targetName', 'seqFrom', 'seqTo', 'Evalue', 'score'], ascending=[True,True, False, False, True])
else:
    melhorResultado = juntoFinal.sort_values(['targetName', 'seqFrom', 'seqTo', 'Evalue', 'score'], ascending=[True,True, False, False, True])

#resultFinal = juntoFinal.sort_values(['targetName', 'seqFrom', 'seqTo', 'Evalue', 'score'], ascending=[True,True, False, False, True])
#resultFinal = resultFinal.reset_index(drop=True)

resultFinal = melhorResultado.reset_index(drop=True)

resultFinal['ncRNA'] = resultFinal['ncRNA'].fillna('ncRNA')
resultFinal['queryName'] = resultFinal['queryName'].fillna('ncRNA')


for i in resultFinal.index:
    resultFinal.loc[i,'ID'] = 'id'+str(i+1)

resultFinalFinal = resultFinal.replace({'minus': '-', 'plus':'+'})
resultFinal.to_csv(nameI+'_annotation.csv', index=False)



gff = salvar_gff(resultFinalFinal, nameI)
qtdTotal_ncRNA(resultFinal, nameI)
graficoTotal_barra(resultFinal, nameI)
graficoTotal_stacked(resultFinal, nameI)
diagrama_venn(resultFinal, nameI)
salvar_original(resultFinal, arqBlast1, arqInfernal)
resultados_tRNA(resultFinal, nameI)
resultados_miRNA(resultFinal, nameI)
if args['fasta'] != None:
    salvar_fasta(gff, nameI)
    

