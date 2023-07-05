######################################################################################################
######################################################################################################
############################################## Implementation of phydrus #############################
######################################################################################################
######################################################################################################
## hydrus model, sampled parameters are taken from the 'spot_setup_hydrus_python.py' module. 
# the results are the modeled soil moisture at the 3 levels (15cm, 40cm and 70cm)
# for details about the parameters refer to (https://github.com/phydrus/phydrus)
 
import os
import pandas as pd
import phydrus as ps
import numpy as np
from distutils.dir_util import copy_tree, remove_tree
import sys


#specify max simulation time and the name of the plot
# tmax=250
# plot="RF36"

def get_moisture(ob_n,tmax,idx):
    #function return the modeled moisture at 15,40 and 70 cm
    #global start,end
    with open(ob_n) as file:
        for i, line in enumerate(file.readlines()):
            if "time" in line:
                start = i-1
            elif "end" in line:
                end = i-1
                break
                
    df = pd.read_csv(ob_n, skiprows=start,nrows=end-start-1,index_col=False,skipinitialspace=True,delim_whitespace=True,error_bad_lines=False)    
    #time_f=idx
    time_f=list(range(1,tmax+1))
    df_f=df[df.time.isin(time_f)]# filter out the data

    return df_f['theta'],df_f['theta.1'],df_f['theta.2']

def check(f):
    #function check if 'end' is spotted in the output
    with open(f) as file:
        if "end" in file.read():
            return True
        return False
def check_nan(h):
    # functioning return true if NAN value spotted in the output file
    with open(h) as file:
        if "NaN" in file.read():
            return True
        return False


def extract_th(g,tmax,idx):
    # first check if simulation was not interrupted by checking if end was written in the output,
    # and second check that there are no NAN values in the output file
    #if fales return negative infinity which will result in high error in the objective function,
    #which lead to searching away from this set of parameters
    if check(g)==True and check_nan(g)==False:
        
        theta15,theta40,theta70 = np.array(get_moisture(g,tmax,idx))
        #theta15=np.array(get_moisture(g).theta)
        #theta40=np.array(get_moisture(g).theta.1)
        #theta70=np.array(get_moisture(g).theta.2)
        ##check if the simulation ran until the end of the end time specified as tmax
        if len(theta15)==tmax:
            return theta15,theta40,theta70 
        
        else:
            theta15=theta40=theta70=np.full((tmax, ), -np.inf)
            return theta15,theta40,theta70
    else:
        theta15=theta40=theta70=np.full((tmax, ), -np.inf)
        return theta15,theta40,theta70
    
def phyd_sim(thetas1,Alfa1,n1,ks1,thetas2,Alfa2,n2,ks2,thetas3,Alfa3,n3,ks3,cpu_id,tmax,idx):
    #copy_tree(self.hymod_path, self.hymod_path+call)
    curdir = os.getcwd()
    #folders copied to a folder with cpu_id on which the model will be run
    copy_tree(curdir, curdir+cpu_id)
    os.chdir(curdir+cpu_id)
    ws = 're_17_test'
    hywd = os.getcwd()
    files=os.listdir(hywd)
    boundary_file=[file for file in files if "boundary" in file][0]
    atmo_file=[file for file in files if "atmo" in file][0]
    # Path to folder containing hydrus.exe 
    #path to hydrus.exe  should be change in Linux workstation
    #hywd="C:\\Users\AlbaraAlmawazreh\\Desktop\\pyhydrus\\phydrus-master\\phydrus-master\\examples\\"
    exe = os.path.join(hywd,'hydrus')  ##needs to be modified
    # Description
    desc = 'RE_13_WATERFLOW_ROOTWATERUPTAKE_ROOTGROWTH'
    # Create model
    ml = ps.Model(exe_name=exe, ws_name=ws, name="model", description=desc,  mass_units="mmol", time_unit="days", length_unit="cm")
    ml.basic_info["lFlux"] = True
    ml.basic_info["lShort"] = False

    ml.add_time_info(tinit=0,tmax=tmax,print_times=False,dt=0.01,dtmin=1*10**-5,dtmax=1)

    #3. Add processes and materials

    ml.add_waterflow(model=3, maxit=10, tolth=0.001, tolh=1, ha=10000, hb=1*10**-6, linitw=True, top_bc=3, bot_bc=4)

    m = ml.get_empty_material_df(n=3)

    #m.loc[[1,2,3]] = [[0.01, 0.31, 0.0295, 1.25, 38, 0.5], [0.0005, 0.36, 0.0246, 1.37, 7, 0.5], [0.00001, 0.395, 0.026, 1.31, 32, 0.5]]
    #m.loc[[1,2,3]] = [[0.01, thetas1, Alfa1, n1, ks1, 0.5], [0.01, thetas2, Alfa2, n2, ks2, 0.5], [0.01, thetas3, Alfa3, n3, ks3, 0.5]]
    m.loc[[1,2,3]] = [[0.01, thetas1, Alfa1, n1, ks1, 0.5], [0.01, thetas2, Alfa2, n2, ks2, 0.5], [0.01, thetas3, Alfa3, n3, ks3, 0.5]]
    
    ml.add_material(m)

    #4. Add profile information¶

    nodes = 450  # Disctretize soil column into n elements
    depth = [-25, -50,-90]  # Depth of the soil column

    ihead = pd.read_csv(boundary_file,decimal=",",sep=";",index_col=0,header=None).index.tolist()

    profile = ps.create_profile(bot=depth, dx=abs(depth[-1] / nodes), h=ihead, mat=m.index)
    profile.iloc[:,1]=ihead
    ml.add_profile(profile)  # Add the profile

    #5. Add observation nodes¶
    ml.add_obs_nodes([-15, -40, -70])

    #6. Add atmosphere boundary conditions¶

    atm = pd.read_csv(atmo_file,decimal=",", sep=";",index_col=0)
    ml.add_atmospheric_bc(atm)

    #7. Add root water uptake¶
    #for labalb
    ml.add_root_uptake(model=0, p0=-10, p2h=-750, p2l=-2000, p3=-8000, r2h=0.5, r2l=0.1, poptm=[-25,-25,-25])
    #for maize(sweet corn)
    #ml.add_root_uptake(model=0, p0=-10, p2h=-500, p2l=-1000, p3=-8000, r2h=0.5, r2l=0.1, poptm=[-25,-25,-25])

    #8. Add root growth¶
    #for lablab plot RF13
    #ml.add_root_growth(irootin=2, irfak=0, trmin=171, trmed=197, trmax=255, xrmin=2, xrmed=11.5, xrmax=90, trperiod=365)
    #for maize plot RF36
    #trmin initial time for root growth, trmed time of known root depth, trmax time of the end of root uptake
    ml.add_root_growth(irootin=2, irfak=0, trmin=69, trmed=103, trmax=148, xrmin=2, xrmed=32.5, xrmax=75, trperiod=365)

    ml.write_input()
    ml.simulate()
    obs=os.path.join(os.getcwd(),ws,"OBS_NODE.OUT")
    th15,th40,th70=extract_th(obs,tmax,idx)
    #filter according to idx in obs
    results=np.concatenate((th15, th40,th70))
    os.chdir(curdir)
    if cpu_id != "opt": remove_tree(curdir+cpu_id)

    return results
    
    
#a=phyd_sim(0.39,1.0595,1.25,50,0.5,0.39,1.0595,1.25,50,0.5,0.39,1.0595,1.25,50,0.5,123,idx)
#ws = "re_17_test" 

#phyd_sim(0.2,0.0295,1.25,38,0.5)
#phyd_sim(0.39,1.0595,1.25,50,0.5,"2132")
#OBS=os.path.join(os.getcwd(),ws,"OBS_NODE.OUT")
#Theta15=np.array(get_moisture(OBS).theta)
#a=Theta15
