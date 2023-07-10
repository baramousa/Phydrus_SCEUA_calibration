# Phydrus_SCEUA_calibration
Calibration of Hydrus-1d model implemented using (Phydrus) using the SCEUA from Spotpy on a linux machine.
The batch file is run which determine the computing resources from a cluster of machines, it compiles the source (hydrus-d1) and then execute the 'excute.py' module. 

The start of the calibration is initiated by running the excute.py module, where one can set the SCEUA parameters, such as the maximum number of iterations and the number of complexes (ngs). Refer to (https://github.com/baramousa/spotpy/blob/master/spotpy/algorithms/sceua.py) for more detailed information.

In phyd.py module the Hydrus-1d implemented using phydrus, is set where it recieve the sampled parameters from the spot_setup_hydrus_python.py module and returns the model output (in this example the soil mositure at certain depths), see (https://github.com/phydrus/phydrus) for more details.

The spot_setup_hydrus_python.py, initiate the calibration drawing the parameters samples and sending them to the phyd.py module, then recieving the results and evaluating them against the observed data. 

Following csv data are needed for running the model calibration:
- atmospheric data (for example 2018_atmo.csv)
- observed soil moisture (such as the RF_observed.csv)
- and boundary conditions, i.e the intial soil moisture in the soil profile (such as the RF_boundary.csv)
