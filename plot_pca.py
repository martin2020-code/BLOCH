import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler
from mpl_toolkits import mplot3d
import os


i = 10
l1, l2, l3 = 1,0.2,0.0025

snapshot = os.listdir('data/VideoData/')[i]
path = os.path.join(
        os.curdir,
        'data', 'VideoData',
        snapshot)
X = pd.read_csv(path, delimiter=' ', header=None).values
sc = StandardScaler(with_std=False)
Xstd = sc.fit_transform(X)


path = os.path.join(
        os.curdir,
        'data', 'PCA',
        snapshot)
C = pd.read_csv(path, delimiter=' ', index_col=None).values[:,1:]


fig = plt.figure()
ax = plt.axes(projection='3d')


ax.plot3D([0,C[0,0]*l1],[0,C[0,1]*l1],[0,C[0,2]*l1], label='pca 1', color='red')
ax.plot3D([0,C[1,0]*l2],[0,C[1,1]*l2],[0,C[1,2]*l2], label='pca 2', color='magenta')
ax.plot3D([0,C[2,0]*l3],[0,C[2,1]*l3],[0,C[2,2]*l3], label='pca 3', color='orange')

ax.scatter(Xstd[:,0],Xstd[:,1],Xstd[:,2])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()
