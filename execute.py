# # -*- coding: utf-8 -*-

# """
#This script initiate the optimisation of hydrus parameters using the SCE-UA algorithm.
#The optimisation are run in parallel, with maximum iterations(simulations) are set to 70000,
# The parameters of SCE-UA were set as follow , number of complexes (ngs=14),kstop=10,peps=0.1, pcento=0.01.
#The meaning of SCE-UA parameters can be found in more details at (https://github.com/baramousa/spotpy/blob/master/spotpy/algorithms/sceua.py) 
# Created on Fri May  7 13:32:42 2021

# @author: AlbaraAlmawazreh
# """


import pandas as pd
import numpy as np
from numpy import inf
import spotpy
import os 
from phyd import phyd_sim
#from spotpy.examples.spot_setup_hymod_exe import spot_setup
from spot_setup_hydrus_python import spot_setup
import matplotlib
import matplotlib.pyplot as plt
from  spotpy.likelihoods import gaussianLikelihoodMeasErrorOut as GausianLike
import math

parallel ='mpi'
spot_setup=spot_setup(spotpy.objectivefunctions.rmse)
sampler=spotpy.algorithms.sceua(spot_setup, dbname='SCEUA_hydrus', dbformat='csv',parallel=parallel)
rep=70000
sampler.sample(rep, ngs=14, kstop=10, peps=0.1, pcento=0.01)  
results=sampler.getdata()

# the following opens the results, calculate percentage of fail model realisation to the total, then sorts all realisation
#from best to worst, and then save the best one percent in a file that can be used to have an idea about the uncertainity of the model
df = pd.DataFrame(data=results, columns=results.dtype.names)
df1=df[df.like1 != inf]
perc_inf=len(df1)/len(df)
print(perc_inf)
len_top=round(len(df1)*1/100)
print(len_top)
len_top=10 if len_top <10 else len_top
df1=df1.sort_values(by=['like1'])
df1=df1.head(len_top)
df1.to_csv('uncertainty.csv',decimal=',',sep=';')
#save df1 to csv
fields=[word for word in df1.columns if word.startswith('sim')]
q5,q95=[],[]
for field in fields:
    q5.append(np.percentile(df1[field],2.5))
    q95.append(np.percentile(df1[field],97.5))

# Figure one plots the iterations vs the error(RMSE)
#figure1
fig= plt.figure(1,figsize=(9,5))
plt.plot(results['like1'])
#plt.show()
plt.ylabel('RMSE')
plt.xlabel('Iteration')
fig.savefig('sce_objectivefunctiontrace.png',dpi=300)


# Figure plot the results of the best model vs the observed soil moisture content
#figure2
tmax=spot_setup.tmax
bestindex,bestobjf = spotpy.analyser.get_minlikeindex(results)
best_model_run = results[bestindex]


best_simulation = list(best_model_run[fields])
#get modeled soil moisture from best_sim for each layer 15,40,70
th15_mod,th40_mod,th70_mod=best_simulation[0:tmax],best_simulation[tmax:2*tmax],best_simulation[2*tmax:3*tmax]
#get the date data
date_DATA= spot_setup.mo
#change the date to datetime
date_DATA.iloc[:,0]=pd.to_datetime(date_DATA.iloc[:,0],format='%Y-%m-%d', errors='ignore')
#filter out the relvent data
#date_DATA=date_DATA[(date_DATA['date'] > '2017-02-18')]

x=date_DATA.iloc[:,0]
x_mod=pd.to_datetime(spot_setup.atmo.date,format='%Y-%m-%d', errors='ignore')
fig2, ax = plt.subplots(3,sharex=True,figsize=(16,9))
#need date observed and date modeled
fig2.subplots_adjust(hspace=0.1)

n_obs=int(len(spot_setup.evaluation())/3)
ax[0].plot(x_mod,th15_mod,color='black',linestyle='solid', label='th15cm')
ax[0].plot(x,spot_setup.evaluation()[0:n_obs],'r.',markersize=3, label='Observed_moisture_15cm')

ax[1].plot(x_mod,th40_mod,color='black',linestyle='solid', label='th40cm')
ax[1].plot(x,spot_setup.evaluation()[n_obs:2*n_obs],'r.',markersize=3, label='Observed_moisture_40cm')

ax[2].plot(x_mod,th70_mod,color='black',linestyle='solid', label='th70cm')
ax[2].plot(x,spot_setup.evaluation()[2*n_obs:3*n_obs],'r.',markersize=3, label='Observed_moisture_70cm')

ax[0].set(title='SOIL MOISTURE BEST EVALUATION VS OBS\nBest objf.='+str(round(bestobjf,3)))
ax[2].set(ylabel="soil moisture [vol.]",
       xlabel="date [days]")
for axes in range(len(ax)):
    ax[axes].legend(loc='upper left')
    ax[axes].set(ylabel="soil moisture [vol.]")

fig2.savefig('sce_best_modelrun.png',dpi=300)
