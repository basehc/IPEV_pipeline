import os
import re
import numpy as np
import subprocess
import time
import sys

f1=open('SRR_list_1.txt','r')

result=[]

for line in f1:
    result.append(line.strip('\n'))


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
for i_name in SRR_fold:
    try:
        dit_name='./'+i_name+'/'
        os.chdir(dit_name)
        if 'before.png' in os.listdir(os.getcwd()):
            print(i_name)
            result.remove(i_name)
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

                    print(store)





        path_parent = os.path.dirname(os.getcwd())
        os.chdir(path_parent)
    except:
        pass
print('error_n is: ',error_n)
print('all_n is: ',all_n)
print('acc is: ',acc)
print(result)
error_srr.sort()
print('error is: ', error_srr)
def main_put():
    wherenow = os.getcwd()
    SRR_fold = file_name(wherenow,'SRR916')
    for i_name in SRR_fold:
        try:
            dit_name='./'+i_name+'/'
            os.chdir(dit_name)
            if len(os.listdir(os.getcwd()))==11:
                subprocess.run(['pkurun-cnlong', '1', '20', 'python', 'cmd_all.py'])
                print(i_name)
            path_parent = os.path.dirname(os.getcwd())
            os.chdir(path_parent)
            

        except:
            print('error')
        time.sleep(0.5)

if str(sys.argv[1])=='1':
    main_put()

def lookup():
    for line in error_srr:
        for check in open('sra_infor.txt'):
            if line in check:
                print(str(check).strip().split()[14])
                
        



if str(sys.argv[2])=='1':
    lookup()