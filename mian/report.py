import os
import re
import numpy as np
import subprocess
import time
import sys

f1=open('SRR_list.txt','r')

result=[]

for line in f1:
    result.append(line.strip('\n'))


f1.close()


def main_put():

    for i_name in result:
        try:
            dit_name='./'+i_name+'/'
            os.chdir(dit_name)

            subprocess.run(['pkurun-cnlong', '1', '20', 'python', 'cmd_all.py'])

            path_parent = os.path.dirname(os.getcwd())
            os.chdir(path_parent)
            

        except:
            print('error')
        time.sleep(0.5)

if str(sys.argv[1])=='1':
    main_put()

