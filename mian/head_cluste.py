import numpy as np 
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import numpy as np
from skbio.stats.distance import anosim
from skbio import DistanceMatrix
from scipy.spatial import distance
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

data=np.load('./all.npy')

pca = PCA(n_components=50)
data=pca.fit_transform(data)
tsne = TSNE(n_components=2,init='pca')


result = tsne.fit_transform(data)
x_min,x_max = np.min(result,0),np.max(result,0)
data = (result-x_min)/(x_max-x_min)
fig = plt.figure()

ax1 = fig.add_subplot(111)

plt.scatter(data[:221,0],data[:221,1],c='blue',s=5,marker='^',label="temperate phage")
plt.scatter(data[221:,0],data[221:,1],c='red',s=5,marker='.',label="virulent phage")
store_infor=dict()

label_x="virulent phage_x"
label_y="virulent phage_y"

store_infor[label_x]=list(data[221:,0])
store_infor[label_y]=list(data[221:,1])
label_x="temperate phage_x"
label_y="temperate phage_y"

store_infor[label_x]=list(data[:221,0])
store_infor[label_y]=list(data[:221,1])
store_infor = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in store_infor.items() ]))
store_infor.to_csv('./before.csv',index=False)

data=distance.cdist(data, data, 'euclidean')
s_label=['s'+str(i) for i in range(data.shape[0])]
dm = DistanceMatrix(data,
                     s_label)
grouping=['Group'+str(i) for i in [1]*220+[2]*int(data.shape[0]-220)]
print(anosim(dm, grouping, permutations=999))



plt.xlabel('tSNE_1')  
plt.ylabel('tSNE_2')

plt.xlim((0, 1))
plt.ylim((0, 1))
my_x_ticks = np.arange(0, 1.1, 0.5)
my_y_ticks = np.arange(0, 1.1, 0.5)
plt.xticks(my_x_ticks)
plt.yticks(my_y_ticks)

fig.savefig('before.png')
