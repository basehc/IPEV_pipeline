import sys
from Bio import SeqIO
look_up = dict()
with open(str(sys.argv[1]),'r') as handle:
    for rec in SeqIO.parse(handle, 'fasta'):
    	look_up[rec.name]=str(rec.seq)
all_label=set(look_up.keys())
#e_value=float(input('input e value:  '))
e_value=0.0001

pro_label=[]
with open(str(sys.argv[2]), 'r') as f1:
    #all_bacteria=[]
    for line in f1.readlines():
        if float(line.split()[10])<e_value:
            #all_bacteria.append(line.split()[1])
            pro_label.append(line.split()[0])

phage_label=all_label-set(pro_label)
phage_name=str(sys.argv[1])+'_phage.fasta'
pro_name=str(sys.argv[1])+'_pro.fasta'

with open(phage_name,'w') as f2:
	for i in phage_label:
		line1='>'+str(i)+'\n'
		f2.writelines(line1)
		line2=look_up[i]+'\n'
		f2.writelines(line2)

with open(pro_name,'w') as f3:
	for i in pro_label:
		line1='>'+str(i)+'\n'
		f3.writelines(line1)
		line2=look_up[i]+'\n'

		f3.writelines(line2)
