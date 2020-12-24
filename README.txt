
============
README
============

This research project develops a solar forecasting model that gives one-day ahead predictions of solar power at the EE building, 
Dept of Electrical Engineering, Faculty of Engineering, Chulalongkorn University

We refer the technical details to

Supachai Suksamosorn, Naebboon Hoonchareon and Jitkomut Songsiri, 
Post-processing of NWP forecasts using Kalman filtering with operational constraints for day-ahead solar power forecasting in Thailand,

http://jitkomut.eng.chula.ac.th/pdf/nwp_moskf_access.pdf

Developers: Supachai Suksamosorn and Jitkomut Songsiri


Input files:

1) wrf_measurement.mat : contain raw measurement data and WRF forecasts 
Data description: EE solar 8kW site from 2017-2018 
2) p_15kw.mat : contain generated power measurement during Apr 2017 - Dec 2018 (15kW site)
3) pv_model.mat : estimated coefficients of PV conversion model (linear model with inputs = Iwrf, Twrf) (8kW solar site)
4) pv_model_15kw.mat : estimated coefficients of PV conversion model (linear model with inputs = Iwrf, Twrf) (15kW solar site)

Main files:

0) graphs_moskf_paper.m
This file plots all the graphs and type tables in the paper.
Required files: simulated results in the folder 'result_files'

To simulate results, run these codes run in order. 
These codes require files in folder 'datainput_files'
1) prepare_model_structure.m : Split data into train and test sets and create model objects
2) run_model.m :  Run our models (MOS, MOS+KF) and models from literature (Pelland, Diagne, Lorenz, Persistent)
3) conversion.m : Calculate predicted power from predicted irradiance from 2). We run this file twice (for two solar sites).
4) misc files (performance_index, performance_index_pv, printtable, printtable_head, rankeer) : analyze the results 

Outputs of the above M files are model objects (structure variables) : 
solar_train, solar_test, MODELNAME_train  and MODELNAME_test
Note: date-time of solar_train, solar_test are time indices of measurement values, 
but date-time of MODELNAME_train, MODELNAME_test are time indices of FORECAST VALUES and residual errors

Example of results in 'result_files' folder 

RESULTS: Forecasting results (Ihat, residual, performance indices) of each model based on validation set approach
daresult_validation1.mat: contain forcasting results based on training from Jan 2017-June 2018 
daresult_validation2.mat: contain forcasting results based on training from Jan 2017-Dec 2017

Moreover, the codes in the folder 'code_all_models' has the hourly model, which is not reported in the paper. 

