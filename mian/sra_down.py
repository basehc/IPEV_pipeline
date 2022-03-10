import wget
import sys
import os
import random
import shutil
from shutil import copyfile
f1=open('SRR_list.txt','r')

result=[]

for line in f1:
    result.append(line.strip('\n'))


f1.close()
for i in result:


    print(i)

    try:
        dit_name='./'+i+'/'
        os.makedirs(dit_name)



        os.chdir(dit_name)   

        sra_id=i 
        url = 'https://sra-pub-run-odp.s3.amazonaws.com/sra/'+ sra_id+'/'+sra_id
        print(url)
        print("Downloading... please wait!!")
        filename = wget.download(url)
        copyfile(filename,str(filename+'.1'))

        os.remove(filename)
        src='../'
        src_files = os.listdir(src)
        for file_name in src_files:
            if file_name[-3:]=='.py':
                full_file_name = os.path.join(src, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name,file_name)
        os.remove('sra_down.py')
        os.remove('clean_memory.py')
        os.remove('report.py')    
        path_parent = os.path.dirname(os.getcwd())
        os.chdir(path_parent)
    except FileExistsError:
    # directory already exists
        print('error')
  
