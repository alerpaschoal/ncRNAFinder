#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 16:08:47 2022

VitorTool version  1.5
Ferramenta V26
RNAcentral version 23
Rfam version 14.10
#Removi o -max_hsps do blast
#FerramentaV28 com funções melhoradas

@author: vitor
"""
import subprocess
import os
import argparse
import datetime as dt
import multiprocessing as mp

print("\nInicio: "+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" \n")

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--fasta", required=True, help= "Fasta file")
ap.add_argument("-o", "--output", required=True, help= "output name")
ap.add_argument("-i", "--pident", type=int, choices=range(0,101), metavar="[0-100]", default= "95", required=False, help= "porcentagem de identidade")
ap.add_argument("-c", "--cover", type=int, choices=range(0,101), metavar="[0-100]", default= "95", required=False, help= "porcentagem de cobertura")
ap.add_argument("-t", "--threads", type=int, default= "1", required=False, help= "number of threads")
ap.add_argument("-r", "--BestHit", required=False, choices=['1', '0'], default= "1", help= "Seleciona o melhor resultado quando apresentar fitas diferentes 1 = Sim ou 0 = Não")
args = vars(ap.parse_args())

path = os.path.realpath(os.path.dirname(__file__))
resultFasta = subprocess.run('UpdateFasta.py -f '+args['fasta'], shell=True).returncode

def execu(cmd):
    subprocess.run(cmd, shell=True)
    
cmd =["blastn -db " +args['output'] +"db -query "+path+"/rnacentral_active.fasta -out blastn_"+args['output']+"_i"+str(args['pident'])+".tsv -num_threads "+str(args['threads'])+" -outfmt '6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send sstrand evalue bitscore qlen slen sframe\ '  -perc_identity "+str(args['pident']), "cmscan -o "+ args["output"]+".out --tblout "+args["output"]+"Infernal.tbl --fmt 2 --cut_ga --rfam --nohmmonly --clanin "+path+"/Rfam.clanin --cpu "+str(args['threads'])+" "+path+"/Rfam.cm " +args["fasta"]]
#cmd =[""]
if resultFasta ==0:
    
    subprocess.run('makeblastdb -in New_'+args['fasta'] +' -input_type fasta -dbtype nucl -title ' +args['output']+ 'db -parse_seqids -out  ' +args['output']+ 'db' , shell=True)
    
    process = []
    for i in cmd:
        p= mp.Process(target=execu, args=(i,))
        process.append(p)
        p.start()
        p.join()
    
    if p.is_alive() ==False:
        subprocess.run("sed '/#/d' "+args["output"]+"Infernal.tbl > aux.tbl", shell=True)
        resultInfernal = subprocess.run("awk '{print $2"+'";"'+"$3"+'";"'+"$4"+'";"'+"$5"+'";"'+"$6"+'";"'+"$7"+'";"'+"$8"+'";"'+"$9"+'";"'+"$10"+'";"'+"$11"+'";"'+"$12"+'";"'+"$13"+'";"'+"$14"+'";"'+"$15"+'";"'+"$16"+'";"'+"$17"+'";"'+"$18"+'";"'+"$19"+'";"'+"$20"+'";"'+"$21"+'";"'+"$22"+'";"'+"$23"+'";"'+"$24"+'";"'+"$25"+'";"'+"$26"+'";"'+"$27}' aux.tbl > "+args["output"]+".tbl", shell=True).returncode
        resultBlast = subprocess.run("awk '{if(($9 - $8 + 1)*100/$15 >= "+str(args['cover'])+") print $0 '} blastn_"+args['output']+"_i"+str(args['pident'])+".tsv > blastn_filtered.tsv", shell=True).returncode

        if resultInfernal ==0 and resultBlast ==0:
            subprocess.run("Pipeline.py -b blastn_filtered.tsv -i "+args["output"]+".tbl -f " + args['fasta'] + " -c "+str(args['cover'])+ " -r "+args['BestHit']+ " -t "+str(args['threads']), shell=True)

            os.remove("FASTA_Temp.pkl")
            os.remove(args["output"]+"db.ndb")
            os.remove(args["output"]+"db.nhr")
            os.remove(args["output"]+"db.nin")
            os.remove(args["output"]+"db.njs")
            os.remove(args["output"]+"db.nog")
            os.remove(args["output"]+"db.nos")
            os.remove(args["output"]+"db.not")
            os.remove(args["output"]+"db.nsq")
            os.remove(args["output"]+"db.ntf")
            os.remove(args["output"]+"db.nto")
            os.remove('New_'+args['fasta'])
            
        os.remove("aux.tbl")

print("\nFim: "+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" \n")


