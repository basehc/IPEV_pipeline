import subprocess
import os

from shutil import copyfile
import shutil



def file_name(file_dir):   
    L=[]  
    for file in os.listdir(file_dir):  
        if os.path.splitext(file)[1] == '.1':  
            L.append(file)      
    return L  
 
wherenow = os.getcwd()
i = file_name(wherenow)[0]


subprocess.run(['./../../sratoolkit.2.11.0-ubuntu64/bin/fastq-dump', '--split-3',i])

file1=i+'_1.fastq'
file2=i+'_2.fastq'
out_file=i+'file'

cmd1_='./../../seqtk/seqtk sample -s100 '+file1+' 400000 > sub1.fq'
cmd2_='./../../seqtk/seqtk sample -s100 '+file2+' 400000 > sub2.fq'


os.system(cmd1_)
os.system(cmd2_)

cmd1_='mv sub1.fq '+file1
cmd2_='mv sub2.fq '+file2

os.system(cmd1_)
os.system(cmd2_)

subprocess.run(["./../../SPAdes-3.15.2-Linux/bin/spades.py", "-t", "30", "-1", file1, '-2', file2, '-o', out_file])
out_file=i+'file'
dir_out=out_file+'/contigs.fasta'
copyfile(dir_out, './contigs.fasta')
subprocess.run(['python','run.py','contigs.fasta'])
subprocess.run(['python','1_split_phage_euk.py'])
subprocess.run(['/appsnew/bioapps/blast-2.9.0+/bin/blastn', '-db', "../../../../all.fasta", "-query", "phage.fasta", "-out", "phage_.fasta.out", "-outfmt", "6", "-task", "blastn", "-max_hsps", "1", "-num_alignments", '1', "-num_threads", "20"])

subprocess.run(['python','create_pro_phage.py',"phage.fasta",'phage_.fasta.out'])
subprocess.run(['python','prophage_weight.py','phage.fasta_pro.fasta'])
subprocess.run(['python','phage_weight.py','phage.fasta_phage.fasta'])

subprocess.run(['python','pro_phage_cluste.py'])
subprocess.run(['python','after_S_anosim.py'])

with open('comb_pha_pro.fasta','wb') as wfd:
    for f in ['phage.fasta_phage.fasta','phage.fasta_pro.fasta']:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)
subprocess.run(['python',"before_weight.py",'comb_pha_pro.fasta'])
subprocess.run(['python',"head_cluste.py"])
