%% Prepare parameters
clear
horizon=10;

load('moskf_train_results'); % predicted irradiance
load('moskf_test_results'); % predicted irradiance
load('solar_splitdata');    % measurement data
load('pv_model_15kw.mat'); % PV model coefficients
load('p_15kw.mat'); % p_15kw.mat contains measurement of solar data (15kW site) during 1 Apr 2017 to Dec 31, 2018

beta=beta15k;

install_cap = 15;

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
% Note (1)
% Since the power data of 15kW site starts from 1 Apr, 2017 
% (while the trained Phat can be predicted from 1 Jan 2017), 
%     we then compute the performance based on data from 1 Apr, 2017
%     This correspond to the array index of 'end-4559' in training set.
%     
% Note (2)
% We neglect the evaluation of the first two days of persistent predictions due to initial
% values in the test data set.
% 
% Note (3)
% We skip the first day of prediction of Diagne, Pelland and our models due to initial value in the test data set

% wrf
[wrf_train.perf_p,wrf_train.perf_p_specific]=performance_index_pv(wrf_train.Phat(end-4559:end),p_train,install_cap);
[wrf_test.perf_p,wrf_test.perf_p_specific]=performance_index_pv(wrf_test.Phat,p_test,install_cap);

% persistence
[persis_train.perf_p,persis_train.perf_p_specific]=performance_index_pv(persis_train.Phat(end-4559:end),p_train,install_cap);
[persis_test.perf_p,persis_test.perf_p_specific]=performance_index_pv(persis_test.Phat(horizon*2+1:end),p_test(horizon*2+1:end),install_cap);

% lorenz
[lorenz_train.perf_p,lorenz_train.perf_p_specific]=performance_index_pv(lorenz_train.Phat(end-4559:end),p_train,install_cap);
[lorenz_test.perf_p,lorenz_test.perf_p_specific]=performance_index_pv(lorenz_test.Phat,p_test,install_cap);

% diagne
[diagne_train.perf_p,diagne_train.perf_p_specific]=performance_index_pv(diagne_train.Phat(end-4559:end),p_train,install_cap);
[diagne_test.perf_p,diagne_test.perf_p_specific]=performance_index_pv(diagne_test.Phat,p_test(horizon+1:end),install_cap);

% pelland
[pelland_train.perf_p,pelland_train.perf_p_specific]=performance_index_pv(pelland_train.Phat(end-4559:end),p_train,install_cap);
[pelland_test.perf_p,pelland_test.perf_p_specific]=performance_index_pv(pelland_test.Phat,p_test(horizon+1:end),install_cap);

% mos daily
[daily_mos_train.perf_p,daily_mos_train.perf_p_specific]=performance_index_pv(daily_mos_train.Phat(end-4559:end),p_train,install_cap);
[daily_mos_test.perf_p,daily_mos_test.perf_p_specific]=performance_index_pv(daily_mos_test.Phat,p_test,install_cap);

% kf1a
[daily_kf1a_train.perf_p,daily_kf1a_train.perf_p_specific]=performance_index_pv(daily_kf1a_train.Phat(end-4559:end),p_train,install_cap);
[daily_kf1a_test.perf_p,daily_kf1a_test.perf_p_specific]=performance_index_pv(daily_kf1a_test.Phat,p_test(horizon+1:end),install_cap);



savefile_train = 'moskf_conversion15kW_train_results';
savefile_test = 'moskf_conversion15kW_test_results';

save(savefile_test,'wrf_test','persis_test','lorenz_test','diagne_test','pelland_test',...
    'daily_mos_test','daily_kf1a_test');

save(savefile_train,'wrf_train','persis_train','lorenz_train','diagne_train','pelland_train',...
    'daily_mos_train','daily_kf1a_train');

%% function for conversion
function Phat = pv_conversion(model,beta)
Ihat=reshape(model.Ihat',[],1);
Twrf=model.Twrf/10;
C=[Ihat Twrf Ihat.*Twrf];
Phat=C*beta;
Phat(Phat<0)=0;
end
function Phat = pv_conversion_wrf(model,beta)
Ihat=model.Iwrf;
Twrf=model.Twrf;
C=[Ihat Twrf Ihat.*Twrf];
Phat=C*beta;
Phat(Phat<0)=0;
end




