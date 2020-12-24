%% README
% This research project develops a solar forecasting model that gives one-day ahead predictions of solar power at the EE building, 
% Dept of Electrical Engineering, Faculty of Engineering, Chulalongkorn University
% 
% We refer the technical details to
% 
% Supachai Suksamosorn, Naebboon Hoonchareon and Jitkomut Songsiri, 
% Post-processing of NWP forecasts using Kalman filtering with operational constraints for day-ahead solar power forecasting in Thailand,
% 
% http://jitkomut.eng.chula.ac.th/pdf/nwp_moskf_access.pdf
% 
% Developers: Supachai Suksamosorn and Jitkomut Songsiri

%% Prepare parameters
clear all;
addpath('./datainput_files/');

% switch between line 6 and 7
SOLAR_PLANT = 1; % 8kW
% SOLAR_PLANT = 2; % 15kW 


if SOLAR_PLANT == 1
    load('pv_model.mat')
    load('moskf_train_results');
    load('moskf_test_results');
    load('solar_splitdata');
    install_cap = 8;
    beta = pv_conv.beta;
    
    savefile_train = 'moskf_conversion8kW_train_results';
    savefile_test = 'moskf_conversion8kW_test_results';
    
elseif SOLAR_PLANT == 2

    conversion15k ; % run this file 
    return
end

horizon=10;

%% PV conversion
% wrf
wrf_train.Phat = pv_conversion_wrf(solar_train,beta);
wrf_test.Phat = pv_conversion_wrf(solar_test,beta);

% persistence
persis_train.Phat = pv_conversion(persis_train,beta);
persis_test.Phat = pv_conversion(persis_test,beta);

% lorenz
lorenz_train.Phat = pv_conversion(lorenz_train,beta);
lorenz_test.Phat = pv_conversion(lorenz_test,beta);

% diagne
diagne_train.Phat = pv_conversion(diagne_train,beta);
diagne_test.Phat = pv_conversion(diagne_test,beta);

% pelland
pelland_train.Phat = pv_conversion(pelland_train,beta);
pelland_test.Phat = pv_conversion(pelland_test,beta);

% mos daily
daily_mos_train.Phat = pv_conversion(daily_mos_train,beta);
daily_mos_test.Phat = pv_conversion(daily_mos_test,beta);

% kf1a
daily_kf1a_train.Phat = pv_conversion(daily_kf1a_train,beta);
daily_kf1a_test.Phat = pv_conversion(daily_kf1a_test,beta);


%% performance index
%     
% Note (1)
% We neglect the evaluation of the two initial days of persistent predictions in the test data set.
% 
% Note (2)
% We neglect the evaluation of the initial date of prediction of Diagne, Pelland and our models in the test data set

% wrf
[wrf_train.perf_p,wrf_train.perf_p_specific]=performance_index_pv(wrf_train.Phat,solar_train.P,install_cap);
[wrf_test.perf_p,wrf_test.perf_p_specific]=performance_index_pv(wrf_test.Phat,solar_test.P,install_cap);

% persistence

[persis_train.perf_p,persis_train.perf_p_specific]=performance_index_pv(persis_train.Phat(horizon*2+1:end),solar_train.P(horizon*2+1:end),install_cap);
[persis_test.perf_p,persis_test.perf_p_specific]=performance_index_pv(persis_test.Phat(horizon*2+1:end),solar_test.P(horizon*2+1:end),install_cap);

% lorenz
[lorenz_train.perf_p,lorenz_train.perf_p_specific]=performance_index_pv(lorenz_train.Phat,solar_train.P,install_cap);
[lorenz_test.perf_p,lorenz_test.perf_p_specific]=performance_index_pv(lorenz_test.Phat,solar_test.P,install_cap);

% diagne
[diagne_train.perf_p,diagne_train.perf_p_specific]=performance_index_pv(diagne_train.Phat,solar_train.P(horizon+1:end),install_cap);
[diagne_test.perf_p,diagne_test.perf_p_specific]=performance_index_pv(diagne_test.Phat,solar_test.P(horizon+1:end),install_cap);

% pelland
[pelland_train.perf_p,pelland_train.perf_p_specific]=performance_index_pv(pelland_train.Phat,solar_train.P(horizon+1:end),install_cap);
[pelland_test.perf_p,pelland_test.perf_p_specific]=performance_index_pv(pelland_test.Phat,solar_test.P(horizon+1:end),install_cap);

% mos daily
[daily_mos_train.perf_p,daily_mos_train.perf_p_specific]=performance_index_pv(daily_mos_train.Phat,solar_train.P,install_cap);
[daily_mos_test.perf_p,daily_mos_test.perf_p_specific]=performance_index_pv(daily_mos_test.Phat,solar_test.P,install_cap);

% kf1a
[daily_kf1a_train.perf_p,daily_kf1a_train.perf_p_specific]=performance_index_pv(daily_kf1a_train.Phat,solar_train.P(horizon+1:end),install_cap);
[daily_kf1a_test.perf_p,daily_kf1a_test.perf_p_specific]=performance_index_pv(daily_kf1a_test.Phat,solar_test.P(horizon+1:end),install_cap);


%% Save data

save(savefile_test,'wrf_test','persis_test','lorenz_test','diagne_test','pelland_test',...
    'daily_mos_test','daily_kf1a_test');

save(savefile_train,'wrf_train','persis_train','lorenz_train','diagne_train','pelland_train',...
    'daily_mos_train','daily_kf1a_train');

%% function for conversion
function Phat = pv_conversion(model,beta)
Ihat=reshape(model.Ihat',[],1);
Twrf=model.Twrf;
C=[Ihat Twrf Ihat.*Twrf];
Phat=C*beta;
Phat(Phat<0)=0;
end
function Phat = pv_conversion_wrf(model,beta)
Ihat=model.Iwrf;
Twrf=model.Twrf*10;
C=[Ihat Twrf Ihat.*Twrf];
Phat=C*beta;
Phat(Phat<0)=0;
end




