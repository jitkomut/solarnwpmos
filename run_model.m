
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

%% User define for run KF 
horizon = 10;
step_hourly = 30*horizon; % for MOS+KF Pelland
tf = 13; % time of forecasting 

%% Predict Ihat from proposed MOS daily-step

% daily model
for i=1:horizon
    daily_mos_train.Ihat(:,i)=daily_mos_train.C(:,:,i)*daily_mos_train.beta(:,i);
    daily_mos_test.Ihat(:,i)=daily_mos_test.C(:,:,i)*daily_mos_test.beta(:,i);
end

%% Predicted solar irradiance from Lorenz model
% Lorenz model does not use current measurement to calculate the prediction
% (No effect from the practical constraint that predictions must be
% released by Tf)

lorenz_train.Ihat = lorenz_train.iwrf - (lorenz_train.C*lorenz_train.beta);
lorenz_test.Ihat = lorenz_test.iwrf - (lorenz_test.C*lorenz_test.beta);

% % Reshape Ihat(T,1) into Ihat(D,h) where T = number of samples, h = forecast horizon, D = number of days 

lorenz_train.Ihat=reshape(lorenz_train.Ihat,horizon,[])';
lorenz_test.Ihat=reshape(lorenz_test.Ihat,horizon,[])';

%% Residual (D,h)

% daily model
daily_mos_test.resid = daily_mos_test.Ihat-reshape(solar_test.I,horizon,[])';

% lorenz
lorenz_train.resid = lorenz_train.Ihat-reshape(lorenz_train.y,horizon,[])';
lorenz_test.resid=lorenz_test.Ihat-reshape(lorenz_test.y,horizon,[])';

% persistence
persis_train.resid = persis_train.Ihat - reshape(solar_train.I,horizon,[])';
persis_test.resid = persis_test.Ihat - reshape(solar_test.I,horizon,[])';

% wrf
wrf_train.resid = wrf_train.Ihat-reshape(solar_train.I,horizon,[])';
wrf_test.resid = wrf_test.Ihat-reshape(solar_test.I,horizon,[])';

%% Performance index
% train
[wrf_train.perf,wrf_train.perf_specific] = performance_index(solar_train.Iwrf,solar_train.I);
[daily_mos_train.perf,daily_mos_train.perf_specific] = performance_index(daily_mos_train.Ihat,solar_train.I);
[lorenz_train.perf,lorenz_train.perf_specific] = performance_index(lorenz_train.Ihat,solar_train.I);

[persis_train.perf,persis_train.perf_specific] = performance_index(persis_train.Ihat(3:end,:),solar_train.I(horizon*2+1:end));
% test
[wrf_test.perf,wrf_test.perf_specific] = performance_index(solar_test.Iwrf,solar_test.I);
[daily_mos_test.perf,daily_mos_test.perf_specific] = performance_index(daily_mos_test.Ihat,solar_test.I);
[lorenz_test.perf,lorenz_test.perf_specific] = performance_index(lorenz_test.Ihat,solar_test.I);

[persis_test.perf,persis_test.perf_specific] = performance_index(persis_test.Ihat(3:end,:),solar_test.I(horizon*2+1:end));

%% Diagne model 
diagne_train = run_kf_hourly(diagne_train,tf);
diagne_train.Ihat = reshape(diagne_train.Iwrf(horizon+1:end),horizon,[])'-diagne_train.Ihat;
diagne_train.Ihat(find(diagne_train.Ihat<0))=0;
diagne_train.time = reshape(diagne_train.time,h,[])' 
diagne_test = run_kf_hourly(diagne_test,tf);
diagne_test.Ihat = reshape(diagne_test.Iwrf(horizon+1:end),horizon,[])'-diagne_test.Ihat;
diagne_test.Ihat(find(diagne_test.Ihat<0))=0;
diagne_test.time = reshape(diagne_test.time,h,[])' 

%% Pelland model 
pelland_train = run_kf_hourly(pelland_train,tf,step_hourly);
pelland_train.Ihat = reshape(pelland_train.Iwrf(horizon+1:end),horizon,[])'-pelland_train.Ihat;
pelland_train.Ihat(find(pelland_train.Ihat<0))=0;
pelland_train.time = reshape(pelland_train.time,h,[])' 
pelland_test = run_kf_hourly(pelland_test,tf,step_hourly);
pelland_test.Ihat = reshape(pelland_test.Iwrf(horizon+1:end),horizon,[])'-pelland_test.Ihat;
pelland_test.Ihat(find(pelland_test.Ihat<0))=0;
pelland_test.time = reshape(pelland_test.time,h,[])' 

%% daily-step MOS+KF1a (Our model reported in the paper)
daily_kf1a_train = run_kf_daily(daily_kf1a_train,tf);
daily_kf1a_train.time=reshape(solar_train.time(h+1:end),m,[])';
daily_kf1a_test = run_kf_daily(daily_kf1a_test,tf);
daily_kf1a_test.time=reshape(solar_test.time(h+1:end),m,[])';
 
%% residual error

% diagne
diagne_train.resid=cal_residual_error(diagne_train.Ihat,solar_train.I,horizon);
diagne_test.resid=cal_residual_error(diagne_test.Ihat,solar_test.I,horizon);
% pelland
pelland_train.resid=cal_residual_error(pelland_train.Ihat,solar_train.I,horizon);
pelland_test.resid=cal_residual_error(pelland_test.Ihat,solar_test.I,horizon);
% kf daily 
daily_kf1a_train.resid=cal_residual_error(daily_kf1a_train.Ihat,solar_train.I,horizon);
daily_kf1a_test.resid=cal_residual_error(daily_kf1a_test.Ihat,solar_test.I,horizon);


%% performance index

% diagne
[diagne_train.perf,diagne_train.perf_specific] = performance_index(diagne_train.Ihat,solar_train.I(horizon+1:end));
[diagne_test.perf,diagne_test.perf_specific] = performance_index(diagne_test.Ihat,solar_test.I(horizon+1:end));
% pelland
[pelland_train.perf,pelland_train.perf_specific] = performance_index(pelland_train.Ihat,solar_train.I(horizon+1:end));
[pelland_test.perf,pelland_test.perf_specific] = performance_index(pelland_test.Ihat,solar_test.I(horizon+1:end));
% kf daily
[daily_kf1a_train.perf,daily_kf1a_train.perf_specific] = performance_index(daily_kf1a_train.Ihat,solar_train.I(horizon+1:end));
[daily_kf1a_test.perf,daily_kf1a_test.perf_specific] = performance_index(daily_kf1a_test.Ihat,solar_test.I(horizon+1:end));


%% Save data

save moskf_test_results wrf_test persis_test daily_mos_test daily_kf1a_test lorenz_test diagne_test pelland_test 
save moskf_train_results wrf_train persis_train daily_mos_train daily_kf1a_train lorenz_train diagne_train pelland_train
save solar_splitdata solar_train solar_test

%% Subfunction: KF hourly model
function kf_hourly=run_kf_hourly(kf_hourly,tf,varargin)
% check variable input arguments
if(length(varargin)>1)
    disp('Error: The input argument is exceed the limit'); return;
    if varargin{1}<=0 
        disp('Error: The step of adaptive noise covariance must be positive'); return;
    end
end
% predefined
A = kf_hourly.A;
C = kf_hourly.C;
W = kf_hourly.W;
V = kf_hourly.V;
zhat_time = kf_hourly.z_init; % zhat(1|0)
P_time = kf_hourly.P_init; % P(1|0)
y = kf_hourly.y;
Ihat = [];
num_sample = size(y,2);
time = kf_hourly.time;
horizon = max(hour(time))-min(hour(time))+1;
vt=[];
wt=[];
hour_left=max(hour(time))-tf;
% run KF to estimate zhat(t+1|t)
for i = 1:num_sample-horizon-hour_left % run from 1 to N-horizon-hour_left iterations and provide predicted value from day 2 until day D 
    yhat(:,i) = C(:,:,i)*zhat_time(:,i);
    % Measurement update
    K(:,:,i) = P_time(:,:,i)*C(:,:,i)'/(C(:,:,i)*P_time(:,:,i)*C(:,:,i)'+V); % kalman gain K(t)=P(t|t-1)C(t)/(C(t)'P(t|t-1)C(t)+V(t)) 
    zhat_mea(:,i) = zhat_time(:,i) + K(:,:,i)*(y(:,i)-yhat(:,i)); % zhat(t|t)
    P_mea(:,:,i) = P_time(:,:,i)-K(:,:,i)*C(:,:,i)*P_time(:,:,i); % P(t|t)

    if ~isempty(varargin)
    % adaptive measurement noise covariance
        if i<=varargin{1}
            vt=[vt y(:,i)-yhat(:,i)];
        else
            vt=vt(:,2:end);
            vt=[vt y(:,i)-yhat(:,i)];
        end
        V=cov(vt'); % 0.006 to 0.04
    % adaptive process noise covariance
        if i<=varargin{1}
            wt=[wt zhat_mea(:,i)-zhat_time(:,i)];
        else
            wt=wt(:,2:end);
            wt=[wt zhat_mea(:,i)-zhat_time(:,i)];
        end
        W=cov(wt'); % 10^-7diag(beta) - 10^-5diag(beta)
    end
        
    % Time update
    zhat_time(:,i+1) = A*zhat_mea(:,i); % zhat(t+1|t)
    P_time(:,:,i+1) = A*P_mea(:,:,i)*A' + W; % P(t+1|t)
    
    % run prediction
    if hour(time(i))==tf   
        for j = 1:horizon
            temp(:,j) = C(:,:,i+hour_left+j)*A^(hour_left+j)*zhat_mea(:,i);   % predicted y(t+h|t) 
        end
        Ihat=[Ihat;temp(1,:)]; % assume the first row of yhat from all model is Ihat
    end
end
% return variables
Ihat(find(Ihat<0))=0; % set yhat<0 = 0
kf_hourly.yhat=yhat;
kf_hourly.Ihat=Ihat;
kf_hourly.zhat_time=zhat_time;
kf_hourly.zhat_mea=zhat_mea;
kf_hourly.P_time=P_time;
kf_hourly.P_mea=P_mea;
kf_hourly.K=K;
end

%% Subfunction: KF daily step
function kf_daily=run_kf_daily(kf_daily,tf,varargin)
% check variable input arguments
if(length(varargin)>1)
    disp('Error: The input argument is exceed the limit'); return;
    if varargin{1}<=0
        disp('Error: The step of adaptive noise covariance must be positive'); return;
    end
end
% predefined
A = kf_daily.A;
C = kf_daily.C;
W = kf_daily.W;
V = kf_daily.V;
horizon = size(C,1);
D = size(C,3); % days = number of day
r = horizon-(max(hour(kf_daily.time(:,1)))-tf); % r = number of hours from the first hour until tf  
F = eye(r,horizon);
zhat_time = kf_daily.z_init;
P_time = kf_daily.P_init;
y = kf_daily.y;
time = kf_daily.time;
zhat_time = reshape(zhat_time,size(zhat_time,1),1,[]);
tf_index = find(hour(time(1:horizon))==tf); % integer index of tf
% run KF to estimate zhat(t+1|t)
for i = 1:D-1 % run from 1 to N-horizon-hour_left iterations and provide predicted value from day 2 until day D-1 
   % Measurement update at time tf
    Kf(:,:,i) = P_time(:,:,i)*C(:,:,i)'*F'/(F*C(:,:,i)*P_time(:,:,i)*C(:,:,i)'*F'+F*V*F'); % kalman gain at time tf Kf(t)=P(t|t-1)C(t)'F'/(FC(t)P(t|t-1)C(t)F'+FV(t)F') 
    zhat_mea_tf(:,i) = zhat_time(:,i) + Kf(:,:,i)*(F*y(:,i)-(F*C(:,:,i)*zhat_time(:,i))); % zhat(t|t) at time tf
    
    % Time update at time tf
    zhat_time_tf(:,i+1) = A(:,:)*zhat_mea_tf(:,i); % zhat(t+1|t) at time tf   
    
    % Forecast at tf for day D+1 
    yhat(:,i+1) = C(:,:,i+1)*zhat_time_tf(:,i+1);
    
    % Measurement update (at the end of day, we have complete measurement,
    % so we follow the regular KF update
    K(:,:,i) = P_time(:,:,i)*C(:,:,i)'/(C(:,:,i)*P_time(:,:,i)*C(:,:,i)'+V); % kalman gain K(t)=P(t|t-1)C(t)'/(C(t)P(t|t-1)C(t)+V(t)) 
    zhat_mea(:,i) = zhat_time(:,i) + K(:,:,i)*(y(:,i)-(C(:,:,i)*zhat_time(:,i))); % zhat(t|t)
    P_mea(:,:,i) = P_time(:,:,i)-K(:,:,i)*C(:,:,i)*P_time(:,:,i); % P(t|t)
    
    % Time update
    zhat_time(:,i+1) = A*zhat_mea(:,i); % zhat(t+1|t)
    P_time(:,:,i+1) = A*P_mea(:,:,i)*A' + W; % P(t+1|t)
end
% return variables
yhat(find(yhat<0))=0; % set yhat<0 = 0
yhat(:,1)=[]; % remove predicted values at day 1
kf_daily.Ihat=yhat';
kf_daily.zhat_time=zhat_time;
kf_daily.zhat_mea=zhat_mea;
kf_daily.zhat_time_tf=zhat_time_tf;
kf_daily.zhat_mea_tf=zhat_mea_tf;
kf_daily.P_time=P_time;
kf_daily.P_mea=P_mea;
kf_daily.K=K;
end

%% Subfunction: residual error
function residual_error = cal_residual_error(predicted,measured,horizon)
% MEASURED IS OF THE FORMAT   'TIME x 1' (column vector)
    residual_error = predicted-reshape(measured(horizon+1:end),horizon,[])';
end



