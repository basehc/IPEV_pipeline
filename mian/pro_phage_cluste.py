import numpy as np
import random
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure()

phage=np.load('./1.npy')
pro=np.load('./2.npy')


print(phage.shape)
print(pro.shape)
long_=phage.shape[0]+pro.shape[0]
data = np.concatenate((pro, phage), axis=0)


tsne = TSNE(n_components=2,init='pca')
result = tsne.fit_transform(data)

num1=220
num2=int(float(pro.shape[0])/float(phage.shape[0])*num1)

phage_array=random.sample(range(phage.shape[0]),num1)

x_min,x_max = np.min(result,0),np.max(result,0)
data = (result-x_min)/(x_max-x_min)

with open('data.npy', 'wb') as f_0:

	np.save(f_0,data)


for i in phage_array:
    plt.scatter(data[i,0],data[i,1],c='blue',s=5,marker='x')

pro_array=random.sample(range(phage.shape[0],long_),num2)
for i in pro_array:
    plt.scatter(data[i,0],data[i,1],c='red',s=5,marker='.')



with open('phage.npy', 'wb') as f_1:

	np.save(f_1,phage_array)


with open('pro.npy', 'wb') as f_2:

	np.save(f_2,pro_array)	



fig.savefig('show.png')