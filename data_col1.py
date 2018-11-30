# coding: utf-8
import wntr
import pandas as pd
import numpy as np
inp_file = 'Net5.inp'
wn = wntr.network.WaterNetworkModel(inp_file)
sim = wntr.sim.EpanetSimulator(wn)
results = sim.run_sim()
pressure_at_5hr = results.node['pressure']
pressure_at = results.node['pressure'].loc[0, :]


pred_list = pd.DataFrame()
pred_list = pred_list.append(pressure_at_5hr)
for junction_name, junction in wn.junctions():    
    junction.demand_timeseries_list[0].base_value = junction.demand_timeseries_list[0].base_value*1.1
    sim = wntr.sim.EpanetSimulator(wn)
    results = sim.run_sim()
    pressure_at_5hr = results.node['pressure']
    pred_list = pred_list.append(pressure_at_5hr)
    junction.demand_timeseries_list[0].base_value = junction.demand_timeseries_list[0].base_value/1.1
#pred_list = pred_list.T
a = pred_list.iloc[0].copy()
for i in range(len(pred_list)):    
    pred_list.iloc[i] = a-pred_list.iloc[i]


#df=df.set_index(['columns'],inplace=True)

pred_list.reset_index(drop = True,inplace=True)
pred_list.drop(pred_list.index[0],inplace=True)
pred_list.drop(['1','2'],axis=1,inplace=True)
datas=[]
for i in range(len(pred_list)):
    data=pred_list.iloc[i,i]
    datas.append(data)
#pred_list=pred_list.append(datas)
pred_list = pred_list.T
#pred_list.loc['0'] = datas
for i in range(len(pred_list)):    
    pred_list.iloc[i] = pred_list.iloc[i]/datas

#data_train=pred_list
#from sklearn.cluster import KMeans
#from sklearn import preprocessing 
#data_train0 = preprocessing.scale(data_train) #标准化
#km = KMeans(n_clusters=6).fit(data_train0)
#data_train['cluster'] = km.labels_
#
#for junction_name, junction in wn.junctions(): 
#    b = junction.coordinates
#    
#a,b,c=[],[],[]
#for junction_name, junction in wn.junctions():    
#    a.append(junction.coordinates)
#
#for i in range(len(a)):
#    b.append(a[i][0])
#    c.append(a[i][1])
#
#import matplotlib.pyplot as plt
#colors = np.array(['red', 'green', 'blue', 'yellow','gold','maroon'])
#plt.scatter(b, c, c=colors[data_train["cluster"]])  
#wntr.graphics.plot_network(wn, node_attribute=pressure_at, node_size=30, 
#                        title='Pressure at 5 hours') 
    