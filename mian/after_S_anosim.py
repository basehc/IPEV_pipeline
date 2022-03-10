import numpy as np
from skbio.stats.distance import anosim
from skbio import DistanceMatrix
from scipy.spatial import distance
phage_array=np.load('./phage.npy')
pro_array=np.load('./pro.npy')
print(pro_array.shape)
print(phage_array.shape)


data=np.load('./data.npy')

'''

grouping = ['Group1', 'Group1', 'Group2', 'Group2']
'''
data = np.concatenate((data[pro_array,:],data[phage_array,:]), axis=0)
data=distance.cdist(data, data, 'euclidean')
s_label=['s'+str(i) for i in range(int(len(pro_array)+len(phage_array)))]
dm = DistanceMatrix(data,
                     s_label)
grouping=['Group'+str(i) for i in [1]*len(pro_array)+[2]*len(phage_array)]
print(anosim(dm, grouping, permutations=999))