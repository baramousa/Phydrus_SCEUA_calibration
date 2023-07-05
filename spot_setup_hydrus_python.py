##################################################################################
#################################################################################
################################ Spotpy sampler#################################
#################################################################################
#################################################################################

# This module is spotpy sampler that draw samples of the parameters from the Normal
# and uniform priors (also known prior distributions) then simulate them and using 
#the objective fuction evaluate them against the observed soil moisture.


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from spotpy.parameter import Uniform
from spotpy.parameter import Normal
from spotpy.objectivefunctions import rmse
import pandas as pd
from phyd import phyd_sim
import os
import sys
import numpy as np
import math

class spot_setup(object):
    
    #thetas  = Normal(mean=0.4 , stddev=0.01)
    #Alfa  = Normal(mean=0.04 , stddev=0.001)
    #n  = Normal(mean=1.3 , stddev=0.01)
    #ks  = Normal(mean=1000 , stddev=1.0)
    #L  = Normal(mean=0.5 , stddev=0.01)
    
    thetas1  = Normal(mean=0.39 , stddev=0.035)
    Alfa1    = Uniform(low=0.0001 , high=0.1)
    n1       = Uniform(low=1.01 , high=1.8)
    ks1      = Uniform(low=0.1 , high=10000)
    #L1       = Uniform(low=0.1 , high=2.5)    
    
    thetas2  = Normal(mean=0.39 , stddev=0.035)
    Alfa2    = Uniform(low=0.0001 , high=0.1)
    n2       = Uniform(low=1.01 , high=1.8)
    ks2      = Uniform(low=0.1 , high=10000)
    #L2       = Uniform(low=0.1 , high=2.5)   
    
    thetas3  = Normal(mean=0.39 , stddev=0.035)
    Alfa3    = Uniform(low=0.0001 , high=0.1)
    n3       = Uniform(low=1.01 , high=1.8)
    ks3      = Uniform(low=0.1 , high=10000)
    #L3       = Uniform(low=0.1 , high=2.5) 
    
    
    def __init__(self,obj_func=None):
        #initialy the observed and atmospheric data are read accordingle the index of the 
        #days of intrest are extracted as well as the tmax(or time at the end of the model)
        self.files=os.listdir(os.getcwd())
        #reading the observed moisture content
        self.obs_file=[file for file in self.files if "observed" in file][0]
        self.atmo_file=[file for file in self.files if "atmo" in file][0]
        self.mo=pd.read_csv(self.obs_file,decimal=",",sep=";")
        self.mo=self.mo.dropna()
        self.atmo=pd.read_csv(self.atmo_file,decimal=",",sep=";")
        self.idx=self.mo.idx.tolist()
        self.tmax=len(self.atmo)
        self.th15_obs=self.mo.iloc[:,1].to_numpy(dtype=float)
        self.th40_obs=self.mo.iloc[:,3].to_numpy(dtype=float)
        self.th70_obs=self.mo.iloc[:,5].to_numpy(dtype=float)
        self.trueObs=np.concatenate((self.th15_obs,self.th40_obs,self.th70_obs),axis=None).tolist()
        
        ##incase temp is compared in the objective function
        # self.t15_obs=self.mo.iloc[:,2].to_numpy(dtype=float)
        # self.t40_obs=self.mo.iloc[:,4].to_numpy(dtype=float)
        # self.t70_obs=self.mo.iloc[:,6].to_numpy(dtype=float)
        # self.trueObs=np.concatenate((self.t15_obs,self.t40_obs,self.t70_obs),axis=None).tolist()
        self.obj_func=obj_func
        
        # #old read observation data and forcing data
        # self.mo = pd.read_csv("theta_obs.csv",decimal=".",sep=";",header=None)
        # self.trueObs=self.mo.iloc[0:765, 1].values.tolist()
        # self.obj_func = obj_func
        
        
    def simulation(self,x):
        # the simulation model extract the cpu id and then start a simulation based on the sampled set of parameters and 
        #model them in the phyd_sim module
        if os.environ.get('OMPI_COMM_WORLD_RANK'):

            cpu_id = str(int(os.environ['OMPI_COMM_WORLD_RANK']))

        elif os.environ.get('SLURM_PROCID'):

            cpu_id = str(int(os.environ['SLURM_PROCID']))
        simulations=[]
        #Here the model is actualy started with a unique parameter combination that it gets from spotpy for each time the model is called
        sim = phyd_sim(x[0], x[1], x[2], x[3], x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],cpu_id,self.tmax,self.idx).tolist()
        simulations.append(sim)

        return simulations[0]

    def evaluation(self):
     
        return self.trueObs

    def objectivefunction(self,simulation,evaluation, params=None):
        #cacluate the objective fuction (RMSE)
        a=simulation
        th15_mod=np.array(a[0:self.tmax])
        th40_mod=np.array(a[self.tmax:2*self.tmax])
        th70_mod=np.array(a[2*self.tmax:3*self.tmax])
        ix=np.array(self.idx)
        ix=ix-1
        sims=np.concatenate((th15_mod[ix], th40_mod[ix],th70_mod[ix])).tolist()
        #SPOTPY expects to get one or multiple values back, 
        #that define the performance of the model run
        #obs_len=int(len(self.trueObs)/3)
        #rmse=math.sqrt((((sims-trueObs)**2).sum())/len(trueObs))
        #rmse15=math.sqrt((((sims[0:obs_len]-self.trueObs[0:obs_len])**2).sum())/len(self.trueObs[0:obs_len]))
        #rmse40=math.sqrt((((sims[obs_len:2*obs_len]-self.trueObs[obs_len:2*obs_len])**2).sum())/len(self.trueObs[obs_len:2*obs_len]))
        #rmse70=math.sqrt((((sims[2*obs_len:3*obs_len]-self.trueObs[2*obs_len:3*obs_len])**2).sum())/len(self.trueObs[2*obs_len:3*obs_len]))
        #PLOTS FROM EXECUTE        
        if not self.obj_func:
            
            # This is used if not overwritten by user
            like = rmse(evaluation,sims)
        else:
            #Way to ensure flexible spot setup class
            like = self.obj_func(evaluation,sims)    
        return like



if __name__=="__main__":
    tst=spot_setup()
    tst.trueObs
    
        

