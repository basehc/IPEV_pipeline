from Bio import SeqIO
from sgt import SGT
import pandas as pd
import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
import numpy as np
import tensorflow as tf
import sys
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 

f=open('log',"r")
lines=f.readlines()
phage=dict()
for x in lines:
    if "Prokaryotes" in x:
        #print(x.split()[1])
        phage[int(float(x.split()[1]))]=x.split()[5]
#print(phage.values())

eukar=dict()

for x in lines:
    if "Eukaryotes" in x:
        #print(x.split()[1])
        eukar[int(float(x.split()[1]))]=x.split()[5]
#print(eukar.values())
print(len(lines),len(eukar.values()),len(phage.values()))
f_1=open('phage.fasta','w')
f_2=open('eukar.fasta','w')

with open('contigs.fasta','r') as handle:
    for seq_id, rec in enumerate(SeqIO.parse(handle, 'fasta')):

        if seq_id in phage.keys():

            f_1.writelines('>'+str(rec.name))
            f_1.writelines('\n')
            f_1.writelines(str(rec.seq))
            f_1.writelines('\n')

        if seq_id in eukar.keys():

            f_2.writelines('>'+str(rec.name))
            f_2.writelines('\n')
            f_2.writelines(str(rec.seq))
            f_2.writelines('\n')
            


