import shutil
import os
#
def file_name(file_dir,str_):   
    L=[]  
    for file in os.listdir(file_dir):  
        if str_ in file:  
            L.append(file)      
    return L  
 
wherenow = os.getcwd()
SRR_fold = file_name(wherenow,'SRR916')
for i in SRR_fold:
    try:
        dit_name='./'+i+'/'
        os.chdir(dit_name)
        if 'before.png' in os.listdir(os.getcwd()):

            file_fold=file_name(os.getcwd(),'file')
            print(file_fold)
            for iii in file_fold:
                shutil.rmtree(iii)
            
            fastq_fold=file_name(os.getcwd(),'fastq')
            print(fastq_fold)
            for iii in fastq_fold:
                os.remove(iii)
        path_parent = os.path.dirname(os.getcwd())
        os.chdir(path_parent)
    except:
        pass
                
