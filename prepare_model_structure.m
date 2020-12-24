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

%% Prepare input data for MOS and MOS+KF both hourly and daily model
% input file is data_17to18.mat and output are solar_train and solar _test 
%% ################ User define ############
clear all; close all; clc;
addpath('./datainput_files/');

start_train = '2017-01-01 00:00:00'; % the begining date and time of training data in the format of yyyy-MM-dd HH:mm:ss
last_train = '2018-06-30 23:59:00'; % the last date and time of training data in the format of yyyy-MM-dd HH:mm:ss
start_test = '2018-07-01 00:00:00'; % the begining date and time of testing data in the format of yyyy-MM-dd HH:mm:ss
last_test = '2018-12-31 18:59:00'; % the last date and time of testing data in the format of yyyy-MM-dd HH:mm:ss
h1 = 7; % the begining hour of the day for training and test
h2 = 16; % the last hour of the day for training and test
tf = 13; % set time of forecast for prepare persistence forecast data
load('wrf_measurement.mat')

%% #########################################
% check out of range of the data
if datenum(start_train)<datenum(measure_downsample.time(1)) || datenum(start_test)<datenum(measure_downsample.time(1)) || datenum(start_train)<datenum(wrf_downsample.time(1)) || datenum(start_test)<datenum(wrf_downsample.time(1))
    disp('Error: Data out of range, please choose the new date time of train and test set that available in the data');
    return
elseif hour(measure_downsample.time(1))~= 0 || hour(wrf_downsample.time(1)) ~= 0
    disp('Error: The first data must begin with 00:00 hrs.');
    return
elseif h2-h1<9
    disp('Error: The proposed model use time interval at least 9 hours');
    return
end

%% index for training set and test set
start_train = datetime(start_train,'InputFormat','yyyy-MM-dd HH:mm:ss');
last_train = datetime(last_train,'InputFormat','yyyy-MM-dd HH:mm:ss');
start_test = datetime(start_test,'InputFormat','yyyy-MM-dd HH:mm:ss');
last_test = datetime(last_test,'InputFormat','yyyy-MM-dd HH:mm:ss');

%% prepare C, y and time index for hourly and daily model
horizon=h2-h1+1;
unit='Irradiance (kW/m^2), Relative humidity (0-1), Temperature (10*C), Solar zenith angle (0-1)';

%% construct training data and test data
solar_train = create_struct_data(start_train,last_train,measure_downsample,wrf_downsample,h1,h2);
solar_test = create_struct_data(start_test,last_test,measure_downsample,wrf_downsample,h1,h2);

%% wrf
wrf_train.Ihat = reshape(solar_train.Iwrf,horizon,[])';
wrf_train.time = reshape(solar_train.time,horizon,[])';
wrf_test.Ihat = reshape(solar_test.Iwrf,horizon,[])';
wrf_test.time=reshape(solar_test.time,horizon,[])';

%% construct persistence model (k_index)
persis_train = create_persistence(start_train,last_train,measure_downsample,h1,h2);
persis_train.k = rearrange_persis([zeros(1,h2-h1+1); reshape(solar_train.I./solar_train.Iclr,horizon,[])'],h1,h2,tf);
persis_train.Ihat = persis_train.k.*reshape(solar_train.Iclr,horizon,[])';
persis_train.Twrf = solar_train.Twrf*10; 
persis_train.time = reshape(solar_train.time,horizon,[])';

persis_test = create_persistence(start_test,last_test,measure_downsample,h1,h2);
persis_test.k = rearrange_persis([zeros(1,h2-h1+1); reshape(solar_test.I./solar_test.Iclr,horizon,[])'],h1,h2,tf);
persis_test.Ihat = persis_test.k.*reshape(solar_test.Iclr,horizon,[])';
persis_test.Twrf = solar_test.Twrf*10; % change the scale and neglect the first two days
persis_test.time = reshape(solar_test.time,horizon,[])';


%% model: MOS Lorenz
% y(t) = bias(t) = Iwrf(t)-I(t) = polynomial order 4 of khat(t) and SZA(t); khat(t) = Iwrf(t)/Iclr(t) is clear sky index, Iclr(t) is solar irradiance under clear sky condition estimated from Ineichen clear sky model
% dimension: y(T,1), ehat(T,1), Iwrf(T,1); T = number of samples
% Note: train/test.time is timestamp of measurement
% training data for MOS
lorenz_train.info.name='MOS_Lorenz';
lorenz_train.info.model='y(t) = bias(t) = Iwrf(t)-I(t) = bias(t) = b1*khat(t)^4 + b2*SZA(t)^4 + b3*khat(t)^3 + b4*SZA(t)^3 + b5*khat(t)^2 + b6*SZA(t)^2+ b7*khat(t) + b8*SZA(t); khat(t) = Iwrf(t)/Iclr(t) is clear sky index',...
                         'Ihat(t) = Iwrf(t)-y(t)';
lorenz_train.info.unit='Irradiance (kW/m^2), Solar zenith angle (0-1), Clear sky index (0-1)';;
lorenz_train.info.description={'C(t) = polynomial order 4 of khat(t) and SZA(t) (no intercept), beta(t) = [b1 b2 b3 b4 b5 b6 b7 b8]^T',...
                          'dimension: y(T,1), C(T,p), z(p,1); T = number of samples, p = number of predictors'};
lorenz_train.khat=solar_train.Iwrf./solar_train.Iclr;
lorenz_train.C=[lorenz_train.khat.^4 solar_train.SZA.^4 lorenz_train.khat.^3 solar_train.SZA.^3 lorenz_train.khat.^2 solar_train.SZA.^2 lorenz_train.khat solar_train.SZA];
lorenz_train.e=solar_train.Iwrf-solar_train.I;
lorenz_train.beta=lorenz_train.C\lorenz_train.e;
lorenz_train.y=solar_train.I;
lorenz_train.Twrf=solar_train.Twrf*10;
lorenz_train.iwrf=solar_train.Iwrf;
lorenz_train.time=reshape(solar_train.time,horizon,[])';

% test data for MOS
lorenz_test.info.name='MOS_Lorenz';
lorenz_test.info.model='y(t) = b1*khat(t)^4 + b2*SZA(t)^4 + b3*khat(t)^3 + b4*SZA(t)^3 + b5*khat(t)^2 + b6*SZA(t)^2+ b7*khat(t) + b8*SZA(t); khat(t) = Iwrf(t)/Iclr(t) is clear sky index';
lorenz_test.info.unit='Irradiance (kW/m^2), Solar zenith angle (0-1), Clear sky index (0-1)';;
lorenz_test.info.description={'C(t) = polynomial order 4 of khat(t) and SZA(t) (no intercept), beta(t) = [b1 b2 b3 b4 b5 b6 b7 b8]^T',...
                          'dimension: y(T,1), C(T,p), z(p,1); T = number of samples, p = number of predictors'};
lorenz_test.khat=solar_test.Iwrf./solar_test.Iclr;
lorenz_test.C=[lorenz_test.khat.^4 solar_test.SZA.^4 lorenz_test.khat.^3 solar_test.SZA.^3 lorenz_test.khat.^2 solar_test.SZA.^2 lorenz_test.khat solar_test.SZA];
lorenz_test.y=solar_test.I;
lorenz_test.Twrf=solar_test.Twrf*10;
lorenz_test.beta=lorenz_train.beta;
lorenz_test.iwrf=solar_test.Iwrf;
lorenz_test.time=reshape(solar_test.time,horizon,[])';

%% model: KF Diagne
% y(t) = bias(t) = Iwrf(t)-I(t) = [1 Iwrf(t) SZA(t)]*[beta1(t) beta2(t) beta3(t)]^T
% A = I
% C(t) = [1 Iwrf(t) SZA(t)]
% z(1|0) = 0
% P(1|0) = 5I
% W = I
% V = 0.01
% dimension: A(p,p), C(1,p,T), W(p,p), V(1,1), y(1,T), z(p,1), P(p,p); p = number of predictors, T = number of samples
% Note: train/test.time is timestamp of measurement
h=horizon; p=3; 
% train data
diagne_train.info.name = 'KF_Diagne';
diagne_train.info.model = 'y(t) = bias(t) = Iwrf(t)-I(t) = [1 Iwrf(t) SZA(t)]*[beta1(t) beta2(t) beta3(t)]^T','Ihat(t) = Iwrf(t) - y(t)';
diagne_train.info.unit = 'Irradiance (kW/m^2), Solar zenith angle (0-1)';
diagne_train.info.description = {'A=I, C(t)=[1 Iwrf(t) SZA(t)], W=I, V=0.01, z(1|0)=0, P(1|0)=5I',...
                            'dimension: A(p,p), C(1,p,T), W(p,p), V(1,1), y(1,T), z(p,1), P(p,p); p = number of predictors, T = number of samples'};
diagne_train.z_init = zeros(p,1);
diagne_train.A = eye(p);
diagne_train.W = eye(p);
diagne_train.V = 0.01;
diagne_train.P_init = eye(p)*5;
diagne_train.C = reshape([ones(length(solar_train.Iwrf),1) solar_train.Iwrf solar_train.SZA]',1,p,[]); 
diagne_train.y = (solar_train.Iwrf-solar_train.I)';
diagne_train.Twrf = solar_train.Twrf(h+1:end)*10;
diagne_train.I = solar_train.I;
diagne_train.Iwrf = solar_train.Iwrf;
diagne_train.time = solar_train.time(h+1:end);  
% test data
diagne_train.info.name = 'KF_Diagne';
diagne_train.info.model = 'y(t) = bias(t) = Iwrf(t)-I(t) = [1 Iwrf(t) SZA(t)]*[beta1(t) beta2(t) beta3(t)]^T','Ihat(t) = Iwrf(t) - y(t)';
diagne_train.info.unit = 'Irradiance (kW/m^2), Solar zenith angle (0-1)';
diagne_train.info.description = {'A=I, C(t)=[1 Iwrf(t) SZA(t)], W=I, V=0.01, z(1|0)=0, P(1|0)=5I',...
                            'dimension: A(p,p), C(1,p,T), W(p,p), V(1,1), y(1,T), z(p,1), P(p,p); p = number of predictors, T = number of samples'};
diagne_test.z_init = diagne_train.z_init;
diagne_test.A = diagne_train.A;
diagne_test.W = diagne_train.W;
diagne_test.V = diagne_train.V;
diagne_test.P_init = diagne_train.P_init;
diagne_test.C = reshape([ones(length(solar_test.Iwrf),1) solar_test.Iwrf solar_test.SZA]',1,p,[]); 
diagne_test.y = (solar_test.Iwrf-solar_test.I)';
diagne_test.Twrf = solar_test.Twrf(h+1:end)*10;
diagne_test.I = solar_test.I;
diagne_test.Iwrf = solar_test.Iwrf;
diagne_test.time = solar_test.time(h+1:end);

%% model: Pelland
% y(t) = bias(t) = Iwrf(t)-I(t) = [1 Iwrf(t)]*[beta1(t) beta2(t)]^T
% A = I
% C(t) = [1 Iwrf(t)]
% z(1|0) = 0
% P(1|0) = 5*10^-5I
% W = 10^-5I
% V = 0.01
% dimension: A(p,p), C(1,p,T), W(p,p), V(1,1), y(1,T), z(p,1), P(p,p); p = number of predictors, T = number of samples
% Note: train/test.time is timestamp of measurement
h = horizon; p = 2; 
% train data
pelland_train.info.name = 'KF_Pelland';
pelland_train.info.model = 'y(t) = [1 Iwrf(t)]*[beta1(t) beta2(t)]^T';
pelland_train.info.unit = 'Irradiance (kW/m^2)';
pelland_train.info.description = {'A=I, C(t)=[1 Iwrf(t)], W=10^-5I, V=0.01, z(1|0)=0, P(1|0)=5*10^-5I',...
                                'dimension: A(p,p), C(1,p,T), W(p,p), V(1,1), y(1,T), z(p,1), P(p,p); p = number of predictors, T = number of samples'};
pelland_train.z_init = zeros(p,1);
pelland_train.A = eye(p);
pelland_train.W = eye(p)*10^-5;
pelland_train.V = 0.01;
pelland_train.P_init = eye(p)*5*10^-5;
pelland_train.C = reshape([ones(length(solar_train.Iwrf),1) solar_train.Iwrf]',1,p,[]); 
pelland_train.y = (solar_train.Iwrf-solar_train.I)';
pelland_train.Twrf = solar_train.Twrf(h+1:end)*10;
pelland_train.I = solar_train.I;
pelland_train.Iwrf = solar_train.Iwrf;
pelland_train.time = solar_train.time(h+1:end);  
% test data
pelland_train.info.name = 'KF_Pelland';
pelland_train.info.model = 'y(t) = [1 Iwrf(t)]*[beta1(t) beta2(t)]^T';
pelland_train.info.unit = 'Irradiance (kW/m^2)';
pelland_train.info.description = {'A=I, C(t)=[1 Iwrf(t)], W=10^-5I, V=0.01, z(1|0)=0, P(1|0)=5*10^-5I',...
                                'dimension: A(p,p), C(1,p,T), W(p,p), V(1,1), y(1,T), z(p,1), P(p,p); p = number of predictors, T = number of samples'};
pelland_test.z_init = pelland_train.z_init;
pelland_test.A = pelland_train.A;
pelland_test.W = pelland_train.W;
pelland_test.V = pelland_train.V;
pelland_test.P_init = pelland_train.P_init;
pelland_test.C = reshape([ones(length(solar_test.Iwrf),1) solar_test.Iwrf]',1,p,[]); 
pelland_test.y = (solar_test.Iwrf-solar_test.I)';
pelland_test.Twrf = solar_test.Twrf(h+1:end)*10;
pelland_test.I = solar_test.I;
pelland_test.Iwrf = solar_test.Iwrf;
pelland_test.time = solar_test.time(h+1:end);

%% daily model 
% model: WRF+MOS
% y_i(d) = b1_i*Iwrf_i(d) + b2_i*RHwrf_i(d) +  b3_i*Twrf_i(d) + b4_i*SZA_i(d) ; d = day and i is number of model corresponding to time between h1 and h2  
% C_i(d) = [Iwrf_i(d) RHwrf_i(d) Twrf_i(d) SZA_i(d)]
% beta_i(d) = [b1_i b2_i b3_i b4_i]^T
% dimension: y(D,h), C(D,p,h), beta(p,h); D = number of days, p = number of predictors, h = forecast horizon
% Note: train/test.time is timestamp of measurement
h = horizon; 
% training data 
daily_mos_train.info.name = 'MOS';
daily_mos_train.info.model = 'y_i(d) = b1_i*Iwrf_i(d) + b2_i*RHwrf_i(d) +  b3_i*Twrf_i(d) + b4_i*SZA_i(d) ; d=day and i is number of model corresponding to time between h1 and h2';
daily_mos_train.info.unit = unit;
daily_mos_train.info.description = {'C_i(d) = [Iwrf_i(d) RHwrf_i(d) Twrf_i(d) SZA_i(d)], beta_i(d) = [b1_i b2_i b3_i b4_i]^T',...
                         'dimension: y(D,h), C(D,p,h), z(p,h); D = number of days, p = number of predictors, h = forecast horizon'};
daily_mos_train.C = reshape_mos_daily([solar_train.Iwrf solar_train.RHwrf solar_train.Twrf solar_train.SZA],h);
daily_mos_train.y = reshape(solar_train.I,h,[])';
daily_mos_train.Twrf = solar_train.Twrf*10;
[daily_mos_train.beta,daily_mos_train.resid,daily_mos_train.covbeta] = estimate_beta_daily(daily_mos_train.C,daily_mos_train.y);
daily_mos_train.time = reshape(solar_train.time,horizon,[])';
% test data 
daily_mos_test.info.name = 'MOS';
daily_mos_test.info.model = 'y_i(d) = b1_i*Iwrf_i(d) + b2_i*RHwrf_i(d) +  b3_i*Twrf_i(d) + b4_i*SZA_i(d) ; d=day and i is number of model corresponding to time between h1 and h2';
daily_mos_test.info.unit = unit;
daily_mos_test.info.description = {'C_i(d) = [Iwrf_i(d) RHwrf_i(d) Twrf_i(d) SZA_i(d)], beta_i(d) = [b1_i b2_i b3_i b4_i]^T',...
                         'dimension: y(D,h), C(D,p,h), z(p,h); D = number of days, p = number of predictors, h = forecast horizon'};
daily_mos_test.C = reshape_mos_daily([solar_test.Iwrf solar_test.RHwrf solar_test.Twrf solar_test.SZA],h);
daily_mos_test.y = reshape(solar_test.I,h,[])';
daily_mos_test.Twrf = solar_test.Twrf*10;
daily_mos_test.beta = daily_mos_train.beta ;
daily_mos_test.time = reshape(solar_test.time,horizon,[])';

%% model: WRF+MOS+KF1a (measurement noise covariance is block diagonal) This is the MAIN method presented in the paper
% [I_1(d) I_2(d) ... I_h(d)]^T = blkdiag([C_1(d) C_2(d) ... C_h(d)])*[z_1(d) z_2(d) ... z_h(d)]^T  ; d=day 
% A = identity, size phxph
% C_i(d) = [Iwrf_i(d) RHwrf_i(d) Twrf_i(d) SZA_i(d)]
% C(d)= blkdiag([C_1(d), C_2(d), ... , C_h(d)])
% beta_i(d) = [b1_i b2_i b3_i b4_i]^T
% z(d) = [beta_1(d) beta_2(d) ... beta_h(d)]^T
% zhat(1|0) = [beta_1LS beta_2LS ... beta_hLS]^T, P(1|0) = diag(Cov(beta_1LS), Cov(beta_2LS), ... , Cov(beta_hLS)), W = 1/1000*diag(zhat(1|0)), V = diag(Cov(v)), v = [I_1(d) - Ihat_mos_1(d), I_2(d) - Ihat_mos_2(d), ... , I_h(d) - Ihat_mos_h(d)] 
% dimension: A(ph,ph), C(h,ph,D), W(ph,ph), V(h,m), y(m,D), z(ph,1), P(ph,ph); h = forecast horizon, m = number of outputs, p = number of predictors, D = number of days
% Note: train/test.time is timestamp of measurement

h=horizon; m=h; p=size(daily_mos_train.C,2); 
% training data 
daily_kf1a_train.info.name = 'MOS+KF1a';
daily_kf1a_train.info.model = '[I_1(d) I_2(d) ... I_h(d)]^T = blkdiag([C_1(d) C_2(d) ... C_h(d)])*[z_1(d) z_2(d) ... z_h(d)]^T  ; d=day ';
daily_kf1a_train.info.unit = unit;
daily_kf1a_train.info.description = {'A = identity matrix, C(d)= blkdiag([C_1(d), C_2(d), ... , C_h(d)]) where C_i(d) = [Iwrf_i(d) RHwrf_i(d) Twrf_i(d) SZA_i(d)], z(d) = [beta_1(d) beta_2(d) ... beta_h(d)]^T where beta_i(d) = [b1_i b2_i b3_i b4_i]^T, zhat(1|0) = [beta_1LS beta_2LS ... beta_hLS]^T, P(1|0) = diag(Cov(beta_1LS), Cov(beta_2LS), ... , Cov(beta_hLS)), W = 1/1000*diag(zhat(1|0)), V = diag(Cov(v)), v = [I_1(d) - Ihat_mos_1(d), I_2(d) - Ihat_mos_2(d), ... , I_h(d) - Ihat_mos_h(d)]',...
                                'dimension: A(ph,ph), C(h,ph,D), W(ph,ph), V(h,m), y(m,D), z(ph,1), P(ph,ph); h = forecast horizon, m = number of outputs, p = number of predictors, D = number of days'};
daily_kf1a_train.z_init=reshape(daily_mos_train.beta,[],1); % beta is p x h  and z_init is ph x 1
daily_kf1a_train.A = eye(p*h);
daily_kf1a_train.W = abs(diag(daily_kf1a_train.z_init)/10000);
daily_kf1a_train.V = diag(diag(cov(daily_mos_train.resid))); % resid is D x h where D is the number of days
daily_kf1a_train.P_init = blkdiag(daily_mos_train.covbeta{1,:}); % covbeta is h-cell arrays of pxp matrices
daily_kf1a_train.C = [solar_train.Iwrf solar_train.RHwrf solar_train.Twrf solar_train.SZA];
T = size(daily_kf1a_train.C,1);
daily_kf1a_train.C = reshape(mat2cell(daily_kf1a_train.C,ones(1,T), p),h,1,[]);
daily_kf1a_train.C = construct_blkdiag_kf1(daily_kf1a_train.C);
daily_kf1a_train.y = reshape(solar_train.I,m,[]);
daily_kf1a_train.Twrf = solar_train.Twrf(h+1:end)*10;
daily_kf1a_train.time = reshape(solar_train.time(h+1:end),m,[]);

% test data 
daily_kf1a_test.info.name = 'MOS+KF1a';
daily_kf1a_test.info.model = '[I_1(d) I_2(d) ... I_h(d)]^T = blkdiag([C_1(d) C_2(d) ... C_h(d)])*[z_1(d) z_2(d) ... z_h(d)]^T  ; d=day ';
daily_kf1a_test.info.unit = unit;
daily_kf1a_test.info.description = {'A = identity matrix, C(d)= blkdiag([C_1(d), C_2(d), ... , C_h(d)]) where C_i(d) = [Iwrf_i(d) RHwrf_i(d) Twrf_i(d) SZA_i(d)], z(d) = [beta_1(d) beta_2(d) ... beta_h(d)]^T where beta_i(d) = [b1_i b2_i b3_i b4_i]^T, zhat(1|0) = [beta_1LS beta_2LS ... beta_hLS]^T, P(1|0) = diag(Cov(beta_1LS), Cov(beta_2LS), ... , Cov(beta_hLS)), W = 1/1000*diag(zhat(1|0)), V = diag(Cov(v)), v = [I_1(d) - Ihat_mos_1(d), I_2(d) - Ihat_mos_2(d), ... , I_h(d) - Ihat_mos_h(d)]',...
                                'dimension: A(ph,ph), C(h,ph,D), W(ph,ph), V(h,m), y(m,D), z(ph,1), P(ph,ph); h = forecast horizon, m = number of outputs, p = number of predictors, D = number of days'};
daily_kf1a_test.z_init = reshape(daily_mos_train.beta,[],1);
daily_kf1a_test.A = daily_kf1a_train.A;
daily_kf1a_test.W = daily_kf1a_train.W;
daily_kf1a_test.V = daily_kf1a_train.V;
daily_kf1a_test.P_init = daily_kf1a_train.P_init;
daily_kf1a_test.C = [solar_test.Iwrf solar_test.RHwrf solar_test.Twrf solar_test.SZA];
D=size(daily_kf1a_test.C,1);
daily_kf1a_test.C = reshape(mat2cell(daily_kf1a_test.C,ones(1,D), repmat(p,1,1)),h,1,[]);
daily_kf1a_test.C = construct_blkdiag_kf1(daily_kf1a_test.C);
daily_kf1a_test.y = reshape(solar_test.I,m,[]);
daily_kf1a_test.Twrf = solar_test.Twrf(h+1:end)*10;
daily_kf1a_test.time = reshape(solar_test.time(h+1:end),m,[]);

%% function to construct the input data
% Reshape the matrix C(T,p) to C(D,p,h)
% where D is the number of day, p is the number of predictors, h is the number of models
function specific_hour=reshape_mos_daily(C,horizon)
    for i=1:horizon
        specific_hour(:,:,i)=C(i:horizon:end,:);
    end
end

% Construct the matrix C which is included the delay time of data 
% e.g. lag 1: the matrix C = [var(t) var(t-1)], lag 2: the matrix C = [var(t) var(t-1) var(t-2)]
function C_with_delay=construct_delay(C,time,lag)
    lag=lag-1;
    start_index=find(datenum(time(1)+1)==datenum(time)); % the data is start with the first hour of day 2 
    C_with_delay=[];
    for i=0:lag
        C_with_delay=[C_with_delay C(start_index-i:end-i,:)];
    end
end

% Construct block diagonal for hourly-step MOS+KF3 (the input data must be in cell format)
function diag_C=construct_blkdiag_kf_hourly(C)
    diag_C=[];
    for i=1:size(C,1)
        temp=blkdiag(C{i,:});
        diag_C(:,:,i)=temp;
    end
end

% Construct block diagonal for daily-step MOS+KF1 (the input data must be in cell format)   
function diag_C=construct_blkdiag_kf1(C)
    diag_C=[];
    for i=1:size(C,3)
        temp=blkdiag(C{:,:,i});
        diag_C(:,:,i)=temp;
    end
end

function re_persis = rearrange_persis(x,h1,h2,tf)
    p1=x(:,1:tf-h1+1);
    p2=x(:,tf-h1+2:end);
    p2=circshift(p2,1);
    p2(1,:)=0;
    re_persis = [p1 p2];
    re_persis(end,:) =[];
end

% solve regression y=X*beta for 'h' submodels
function [beta,resid,covbeta]=estimate_beta_daily(X,y)
% beta(p,h), resid(T,h), covbeta is a cell array of pxp matrices
% where p is the number of predictors and h is the number of regression models
% and T is the number of samples used in the regression
[T,p,horizon] = size(X);
    for i=1:horizon
        beta(:,i)=X(:,:,i)\y(:,i);
        resid(:,i) = y(:,i)-X(:,:,i)*beta(:,i);
        covbeta{i}=(1/(T-p))*norm(resid(:,i))^2*((X(:,:,i)'*X(:,:,i))\eye(p));
    end
end

%% function for determine train and test data
% find index train and test data
function index=find_index(start_train,last_train,datetime_vector,h1,h2)
    index = find(datenum(start_train)<=datenum(datetime_vector)&datenum(datetime_vector)<=datenum(last_train));
    index_hour = find(h1<=hour(datetime_vector)&hour(datetime_vector)<=h2);
    index = intersect(index,index_hour);
end

% create structure data 
function struct_data = create_struct_data(start_date,last_date,measure_downsample,wrf_downsample,h1,h2)

index_measure = find_index(start_date,last_date,measure_downsample.time,h1,h2);
index_forecast = find_index(start_date,last_date,wrf_downsample.time,h1,h2);
struct_data.info.variable= {'Solar irradiance (I)','Solar power (P)','Solar irradiance from WRF with spatial averaging 49 grid points (Iwrf)',...
    'Relative humidity from WRF with spatial averaging 49 grid points (RHwrf)',...
                           'Temperature from WRF with spatial averaging 49 grid points (Twrf)','Solar zenith angle from Ineichen model (SZA)'};
struct_data.info.unit = 'Irradiance (kW/m^2) Power (kW) Relative humidity (0-1) Temperature (10*C) Solar zenith angle (0-1)';
struct_data.info.nighttime = 'removed';
struct_data.info.sampling = strcat('sampling period of 60 minutes',mat2str(h1),'to',mat2str(h2),'hrs.',datestr(start_date,' dd-mmm-yyyy'),' to',datestr(last_date,' dd-mmm-yyyy'));
struct_data.time = measure_downsample.time(index_measure);
struct_data.I = measure_downsample.I(index_measure)/1000;
struct_data.P = measure_downsample.P(index_measure);
struct_data.Iwrf = wrf_downsample.Iwrf(index_forecast)/1000;
struct_data.RHwrf = wrf_downsample.RHwrf(index_forecast)/100;
struct_data.Twrf = wrf_downsample.Twrf(index_forecast)/10;
struct_data.SZA = cal_sza(struct_data.time);
struct_data.SZA(find(struct_data.SZA<0))=0;
struct_data.Iclr = cal_iclr(struct_data.SZA)/1000;
end

%% create persistence 
function struct_data = create_persistence(start_date,last_date,measure_downsample,h1,h2)
index_measure = find_index(start_date,last_date,measure_downsample.time,h1,h2);
struct_data.info.variable = {'Solar irradiance (I)','Solar power (P)'};
struct_data.info.unit = 'Irradiance (kW/m^2) and Power (kW)';
struct_data.info.nighttime = 'removed';
struct_data.info.sampling = strcat('sampling period of 60 minutes',mat2str(h1),'to',mat2str(h2),'hrs.',datestr(start_date,' dd-mmm-yyyy'),' to',datestr(last_date,' dd-mmm-yyyy'));
struct_data.time = measure_downsample.time(index_measure);
struct_data.I = measure_downsample.I(index_measure)/1000;
struct_data.P = measure_downsample.P(index_measure);
end

%% solar zenith angle calculation
function sza=cal_sza(datetime_vector)
%% ############# prepare parameters ###############
phi = 13.737;              %<<<<<<<<<<< input latitude of EECU site 
lon = 100.532;             %<<<<<<<<<<< input longitude of EECU site 
%% cal solar zenith angle (sza)
doy = day(datetime_vector,'dayofyear');
theta = 2*pi*(doy-1)/365;
delta = (0.006918-0.399912*cos(theta)+0.070257*sin(theta)-0.006759*cos(2*theta) ... 
    +0.000907*sin(2*theta)-0.002697*cos(3*theta)+0.00148*sin(3*theta))*(180/pi); % solar declination in radian (angle between the sun and equator at noon) 
time_vector = hour(datetime_vector)+minute(datetime_vector)/60; %cal solar hour
et = 9.87*sin(2*2*pi*(doy-81)/364)-7.53*cos(2*pi*(doy-81)/364)-1.5*sin(2*pi*(doy-81)/364); % cal equation of time
omega = abs(time_vector+(et/60)-((105-lon)/15)-12)*15; % angular of solar time in degree (15 degree = 1 hour)
sza = cosd(phi)*cosd(delta).*cosd(omega)+sind(phi)*sind(delta); % solar zenith angle (sza)= cos\theta(t)
end

%% solar irradiance under clear sky condition (Iclr) using Ineichen clear sky model
function iclr=cal_iclr(sza)
%% ############# prepare parameters ###############
Isc = 1366.1; % solar constant
alt = 35; % altitude above mean sea level  
TL1 = 4.86; % Linke turbidity factor obtained from Least-squares (customized to Thailand)
f1 = exp(-alt/8000); % parameter of Ineichen clear sky model 
f2 = exp(-alt/1250); % parameter of Ineichen clear sky model 
a1 = (alt*5.09e-5)+0.868; % parameter of Ineichen clear sky model 
a2 = (alt*3.92e-5)+0.0387; % parameter of Ineichen clear sky model 
AM = 1./(sza+0.15*(93.885-acosd(sza)).^-1.253); % air mass calculation

% calculate clear-sky irradiance 
iclr = a1*Isc*sza.*exp(-a2*AM*(f1+f2*(TL1-1)));
end