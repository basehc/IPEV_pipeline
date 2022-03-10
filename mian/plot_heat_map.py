import sys
import os
import random

import numpy as np
import seaborn as sns
sns.set_theme()
import numpy as np
import random
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


########################################
import os
import re
import numpy as np
import subprocess
import time
import sys

f1=open('SRR_list_1.txt','r')


f1.close()


def file_name(file_dir,str_):   
    L=[]  
    for file in os.listdir(file_dir):  
        if str_ in file:  
            L.append(file)      
    return L  
 
wherenow = os.getcwd()
SRR_fold = file_name(wherenow,'SRR916')
error_n=0
all_n=0
acc=0
error_srr=[]
key_store_1 = dict()

key_store_2 = dict()
for i_name in SRR_fold:
    try:
        dit_name='./'+i_name+'/'
        os.chdir(dit_name)
        if 'before.png' in os.listdir(os.getcwd()):
            print(i_name)

            phage=np.load('./1.npy')
            pro=np.load('./2.npy')
            num1=220
            num2=int(float(pro.shape[0])/float(phage.shape[0])*num1)

            print(num1,num2,float(pro.shape[0])/float(phage.shape[0]))

            file_fold=file_name(os.getcwd(),'pyt')
            for iii in file_fold:

                if '.out' in iii:
                    store=[]
                    for line in open(iii):
                        
                        if "p-value" in line:
                            print(line.strip())
                            store.append(float(line.strip().split()[1]))
                    #close(iii)


                    print(store)
                    key_store_1[i_name]=store





        path_parent = os.path.dirname(os.getcwd())
        os.chdir(path_parent)
    except:
        pass
T8_4000_sra=[]
def lookup():
    wherenow = os.getcwd()
    SRR_fold = file_name(wherenow,'SRR916')

    for line in SRR_fold:
        for check in open('sra_infor.txt'):
            if line in check and 'Illumina HiSeq 2500' not in check:
                T8_4000_sra.append(line)

                key_store_2[tuple(str(check).strip().split()[14].split(',')[3].split('_')[2:])]=key_store_1[line]

            



lookup()

########################################

all_pre_store = sorted(key_store_2.keys(), key=lambda t: t[0])



print('ss',key_store_2)
for i in sorted(key_store_2.keys(), key=lambda t: t[0]):
    print(i,key_store_2[i])
print(T8_4000_sra)

new_sra=set(SRR_fold)-set(T8_4000_sra)

new_key_store_1=dict()
for i_name in list(new_sra):
    try:
        dit_name='./'+i_name+'/'
        os.chdir(dit_name)
        if 'before.png' in os.listdir(os.getcwd()):
            print(i_name)

            phage=np.load('./1.npy')
            pro=np.load('./2.npy')
            num1=220
            num2=int(float(pro.shape[0])/float(phage.shape[0])*num1)

            print(num1,num2,float(pro.shape[0])/float(phage.shape[0]))

            file_fold=file_name(os.getcwd(),'pyt')
            for iii in file_fold:

                if '.out' in iii:
                    all_n+=1
                    store=[]
                    for line in open(iii):
                        
                        if "p-value" in line:
                            print(line.strip())
                            store.append(float(line.strip().split()[1]))
                    #close(iii)
                    if store[0]>store[1]:
                        error_n+=1
                        error_srr.append(i_name)
                    if store[1]>=0.05 :
                        acc+=1

                    new_key_store_1[i_name]=store





        path_parent = os.path.dirname(os.getcwd())
        os.chdir(path_parent)
    except:
        pass
print('error_n is: ',error_n)
print('all_n is: ',all_n)
print('acc is: ',acc)
error_srr.sort()
print('error is: ', error_srr)

print(new_key_store_1)
new_key_store_2=dict()
def lookup1():


    for line in new_sra:
        for check in open('sra_infor.txt'):
            if line in check :
                T8_4000_sra.append(line)

                new_key_store_2[tuple(str(check).strip().split()[14].split(',')[3].split('_')[2:])]=new_key_store_1[line]

            



lookup1()

new_key_store_2_all = sorted(new_key_store_2.keys(), key=lambda t: t[0])

for i in range(916,926):
    new=[]
    for j in new_key_store_2_all:
        if str(i)==j[0]:

            new.append(j)
    new.sort(key=lambda x:int(x[1][1:]))
    for ii in new:
        print(new_key_store_2[ii][1])
    print('\n')

for i in range(916,926):
    new=[]
    for j in new_key_store_2_all:
        if str(i)==j[0]:

            new.append(j)
    new.sort(key=lambda x:int(x[1][1:]))
    for ii in new:
        print(ii)
    print('\n')
