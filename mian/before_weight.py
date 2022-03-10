from Bio import SeqIO
from sgt import SGT
import pandas as pd
import os

import numpy as np

import sys
import time
import matplotlib
import pickle
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
start_time = time.time()
t = time.localtime()
current_time = time.strftime("%H_%M_%S_", t)
pro_num=np.load("./pro.npy")
phage_num=np.load("./phage.npy")
new_all=np.concatenate((phage_num,pro_num), axis=None)

file_0 = str(current_time)+str(sys.argv[1])+'_temper_file/'

os.makedirs(file_0)

def getKmers(sequence, size):
    return [sequence[x:x+size] for x in range(len(sequence) - size + 1)]

def sgtfunc(filename,merlen=3):
    listapl=["AAA", "AAC", "AAG", "AAT",
"ACA", "ACC", "ACG", "ACT",
"AGA", "AGC", "AGG", "AGT",
"ATA", "ATC", "ATG", "ATT",
"CAA", "CAC", "CAG", "CAT",
"CCA", "CCC", "CCG", "CCT",
"CGA", "CGC", "CGG", "CGT",
"CTA", "CTC", "CTG", "CTT",
"GAA", "GAC", "GAG", "GAT",
"GCA", "GCC", "GCG", "GCT",
"GGA", "GGC", "GGG", "GGT",
"GTA", "GTC", "GTG", "GTT",
"TAA", "TAC", "TAG", "TAT",
"TCA", "TCC", "TCG", "TCT",
"TGA", "TGC", "TGG", "TGT",
"TTA", "TTC", "TTG", "TTT"]
    dictuncertain = {
        'N':'G',
        'X':'A',
        'H':'T',
        'M':'C',
        'K':'G',
        'D':'A',
        'R':'G',
        'Y':'T',
        'S':'C',
        'W':'A',
        'B':'C',
        'V':'G',
    }

    num = 0


    with open(filename,'r') as handle:


        key_store = dict()
        count_1=0
        count_2=0
        count_3=0
        count_4=0
        
        for seq_id, rec in enumerate(SeqIO.parse(handle, 'fasta')):
            address1=str(rec.seq)



            for word1, word2 in dictuncertain.items():



                address1 = address1.replace(word1, word2)


            sequence_id = np.zeros(5)
            sequence_id[0]=seq_id
            _length_sequence = len(address1)

            if 2>len(address1[_length_sequence//1800*1800:]):
                sequence_id[4]=_length_sequence//1800
            else:
                sequence_id[4]=_length_sequence//1800+1
            for _num_i in range(_length_sequence//1800+1):
                
                sequence_id[1]=_num_i
                if len(address1[_num_i*1800 : _num_i*1800+1800])<2:
                    pass

                elif 2<=len(address1[_num_i*1800 : _num_i*1800+1800])<400:
                    count_1+=1
                    sequence_id[2]=len(address1[_num_i*1800 : _num_i*1800+1800])
                    sequence_id[3]=1

                    key_store[sequence_id.tobytes()]=[count_1, getKmers(address1[_num_i*1800 : _num_i*1800+1800],merlen)]
                elif 400<=len(address1[_num_i*1800 : _num_i*1800+1800])<800:
                    count_2+=1
                    sequence_id[2]=len(address1[_num_i*1800 : _num_i*1800+1800])
                    sequence_id[3]=2
                    key_store[sequence_id.tobytes()]=[count_2, getKmers(address1[_num_i*1800 : _num_i*1800+1800],merlen)]
                elif 800<=len(address1[_num_i*1800 : _num_i*1800+1800])<1200:
                    count_3+=1
                    sequence_id[2]=len(address1[_num_i*1800 : _num_i*1800+1800])
                    sequence_id[3]=3
                    key_store[sequence_id.tobytes()]=[count_3, getKmers(address1[_num_i*1800 : _num_i*1800+1800],merlen)]
                elif 1200<=len(address1[_num_i*1800 : _num_i*1800+1800])<=1800:
                    count_4+=1
                    sequence_id[2]=len(address1[_num_i*1800 : _num_i*1800+1800])
                    sequence_id[3]=4
                    key_store[sequence_id.tobytes()]=[count_4, getKmers(address1[_num_i*1800 : _num_i*1800+1800],merlen)]
        
        global key_store_1 
        global key_store_2 
        global key_store_3 
        global key_store_4 
        global data_set_1
        global data_set_2
        global data_set_3
        global data_set_4
        key_store_1 = dict()
        key_store_2 = dict()
        key_store_3 = dict()
        key_store_4 = dict()
        for key_id in key_store.keys():

            if np.frombuffer(key_id)[3]==1:
                key_store_1[key_id]=key_store[key_id]
            if np.frombuffer(key_id)[3]==2:
                key_store_2[key_id]=key_store[key_id]
            if np.frombuffer(key_id)[3]==3:
                key_store_3[key_id]=key_store[key_id]
            if np.frombuffer(key_id)[3]==4:
                key_store_4[key_id]=key_store[key_id]


    sgt = SGT(alphabets=listapl,kappa=1, 
              flatten=True, 
              lengthsensitive=False, 
              mode='multiprocessing')            




    corpus = pd.DataFrame(key_store_1.values(),columns=['id', 'sequence'])
    sgtpd = sgt.fit_transform(corpus)
    sgtpd.pop('id')
    data_set_1 = sgtpd.values

    corpus = pd.DataFrame(key_store_2.values(),columns=['id', 'sequence'])
    sgtpd = sgt.fit_transform(corpus)

    sgtpd.pop('id')
    data_set_2 = sgtpd.values

    corpus = pd.DataFrame(key_store_3.values(),columns=['id', 'sequence'])
    sgtpd = sgt.fit_transform(corpus)

    sgtpd.pop('id')
    data_set_3 = sgtpd.values

    corpus = pd.DataFrame(key_store_4.values(),columns=['id', 'sequence'])
    sgtpd = sgt.fit_transform(corpus)

    sgtpd.pop('id')
    data_set_4 = sgtpd.values

    return (count_1, count_2, count_3, count_4)


##########################################################################

##########################################################################

####################################################################

(count_1, count_2, count_3, count_4)=sgtfunc(str(sys.argv[1]))




b1=data_set_1
b2=data_set_2
b3=data_set_3
b4=data_set_4

pre_store_1 = dict()
for num_i, key_id in enumerate(key_store_1.keys()):
    pre_store_1[key_id] = b1[num_i]

pre_store_2 = dict()
for num_i, key_id in enumerate(key_store_2.keys()):
    pre_store_2[key_id] = b2[num_i]

pre_store_3 = dict()
for num_i, key_id in enumerate(key_store_3.keys()):
    pre_store_3[key_id] = b3[num_i]

pre_store_4 = dict()
for num_i, key_id in enumerate(key_store_4.keys()):
    pre_store_4[key_id] = b4[num_i]

all_pre_store = {**pre_store_1,**pre_store_2,**pre_store_3,**pre_store_4}






sort_dict_pre = sorted(all_pre_store.keys(), key=lambda t: np.frombuffer(t)[0])

all_count = count_1+count_2+count_3+count_4
start = 0 


n=0
head=np.zeros([1,64*64])
while start<all_count:
    n+=1
    step = int(np.frombuffer(sort_dict_pre[start])[4])
    denominator = 0
    weight = 0
    for i in sort_dict_pre[start:start+step]:
        denominator += np.frombuffer(i)[2]
        weight += np.frombuffer(i)[2]*all_pre_store[i]
    score = weight/denominator
    if n in new_all:
        print(n)
        score=np.expand_dims(score,axis=0)
        head=np.append(head,score, axis=0)
    
    start +=step

time_run = time.time() - start_time
print("PEVI run time --- %s seconds ---" % (time_run))


head=np.delete(head,0,axis=0)
with open('all.npy', 'wb') as f:
    np.save(f, head)


