# coding: utf-8
#import pyodbc  
##from data_smooth1 import savgol
import numpy as np
import pandas as pd
#import time
#import datetime
#import random
#import tensorflow as tf
#from scipy import interpolate 


data0 = pd.read_excel('邻接矩阵.xlsx',sheet_name= 0) #pd.read_excel默认生成DataFrame对象
data1 = pd.read_excel('邻接矩阵.xlsx',sheet_name= 1) #pd.read_excel默认生成DataFrame对象
data2 = pd.read_excel('邻接矩阵.xlsx',sheet_name= 2) #pd.read_excel默认生成DataFrame对象
data3 = pd.read_excel('邻接矩阵.xlsx',sheet_name= 3) #pd.read_excel默认生成DataFrame对象

#data0.loc[:,:] = 0
index = data1['ID'].tolist()
column = data2['ID'].tolist()
zero_matrix0  = np.random.randint(0,1,(91,112))
data00 = pd.DataFrame(data=zero_matrix0,index=index,columns=column)
##data00[0]
##print(data00[20][0])
##print(data2['ID'])

for i in range(len(data2)):
    data00[data2['ID'][i]][data2['Node1'][i]] = 1
    data00[data2['ID'][i]][data2['Node2'][i]] = -1
A = data00[0:-2]
data01 = 1000*data3['长度']/(data3['流量']*data3['单位水头损失']*1.852)
data_pump = 1/(2*1606.26*7.251e-5)
zero_matrix1 = np.random.randint(0,1,(112,112))
data02 = pd.DataFrame(data=zero_matrix1,index=column,columns=column)

data03 = data01.reset_index()

for i in range(len(data03)):
    data02[data03['index'][i]][data03['index'][i]] = data03[0][i]

data02[3][3] = data_pump

a = np.mat(A)
b = np.mat(data02)
c=-(a*b*(a.T)).I
data_train = pd.DataFrame(data=c,index=index[0:-2],columns=index[0:-2])

from sklearn.cluster import KMeans
from sklearn import preprocessing 
data_train0 = preprocessing.scale(data_train) #标准化
km = KMeans(n_clusters=4).fit(data_train0)
data_train['cluster'] = km.labels_
data_show = data_train['cluster']
data0['cluster'] = data_show
centers = data0.groupby("cluster").mean().reset_index()

import matplotlib.pyplot as plt
colors = np.array(['red', 'green', 'blue', 'yellow'])
plt.scatter(data0["x"], data0["y"],c=colors[data0["cluster"]])

#plt.scatter(centers.x, centers.y, linewidths=3, marker='+', s=300, c='black')
#data3['长度'][335] = 1000
#'red', 'green', 'blue', 'yellow','black','gold','orange','maroon','lightgray'
#data01=data2['flow']
#data02=data2['date']
#data03=data2.index
#'''保存为.xls文件'''
#dataframe = pd.DataFrame({'date':ylist,'flow':datas})    #列表转化成dataframe
#dataframe.to_excel('dataQ8.xlsx',index=False)































#a=np.mat([[1,1,0,0,-1,0],
#          [0,-1,1,0,0,0],
#          [0,0,0,-1,1,0],
#          [0,0,-1,1,0,1]])
#b=np.mat([[3.1094,0,0,0,0,0],
#          [0,10.5785,0,0,0,0],
#          [0,0,6.9982,0,0,0],
#          [0,0,0,8.0216,0,0],
#          [0,0,0,0,16.7184,0],
#          [0,0,0,0,0,2.2667]])
#
#c=-(a*b*(a.T)).I
#print(c)
#print(b)










#import networkx as nx #导入NetworkX包,为了少打几个字母,将其重命名为nx 
#G = nx.Graph() #建立一个空的无向图G 
#G.add_node(1) #添加一个节点1 
#G.add_edge(2,3) #添加一条边2-3(隐含着添加了两个节点2、3) 
#G.add_edge(3,2) #对于无向图,边3-2与边2-3被认为是一条边 
#a = G.nodes()
#print (a) #输出全部的节点: [1, 2, 3] 
#print (G.edges()) #输出全部的边:[(2, 3)] 
#print (G.number_of_edges()) #输出边的数量:1
#print (G.edges(3))
#H = G.to_undirected()
#print (H.nodes()) #输出全部的节点: [1, 2, 3] 
#print (H.edges()) #输出全部的边:[(2, 3)] 
#print (H.number_of_edges()) #输出边的数量:1
#print (H.edges(3))
#G.add_weighted_edges_from([(0,1,3.0),(1,2,7.5)]) 
#print (G.get_edge_data(1,2)) #输出{'weight': 7.5}
#print (G.edges()) #输出全部的边:[(2, 3)] 
#path=nx.all_pairs_shortest_path(G) #调用多源最短路径算法,计算图G所有节点间的最短路径 
#print (path[0][2]) #输出节点0、2之间的最短路径序列: [0, 1, 2] 


#G = nx.Graph()
#nx.add_path(G, [0, 1, 2])
#nx.add_path(G, [0, 10, 2])
#print([p for p in nx.all_shortest_paths(G, source=0, target=2)])