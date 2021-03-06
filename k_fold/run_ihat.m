clear all;

addpath('./input_data/');
% VARIABLE LIST IN 'KFOLD, k_fold_test' = 1) I 2) Iwrf 3) RHwrf 4) Twrf 5) SZA 6) T 7) Iclr 8) P

load('k_fold_data17to18.mat')

savefilepath = './result_files/kfold_ihat.mat' ;


cv = cvpartition(kfold_date,'KFold',10,'Stratify',false);
split_data = mat2cell(kfold,cv.TestSize);
split_date = mat2cell(kfold_date,cv.TestSize);
num_fold =length(split_data);
index=[1:num_fold];

%% prepare data
for i=1:num_fold
    merge_index = setdiff(index,i);
    k_fold_train(:,:,i) = cell2mat(split_data(merge_index));
    k_fold_test(:,:,i) = cell2mat(split_data(i));
    k_fold_train_date(:,i) = split_date(merge_index);
    k_fold_test_date(:,i) = split_date(i);
end
k_fold_train_date=[k_fold_train_date{:}];
k_fold_test_date=[k_fold_test_date{:}];

measured_I = squeeze(k_fold_test(:,1,:)); % 730 x num_fold
ihat_wrf = squeeze(k_fold_test(:,2,:)); % 730 x num_fold
temp_wrf = 10*squeeze(k_fold_test(:,4,:)); % read the temperature and change the unit
p8k = squeeze(k_fold_test(:,8,:)); 

%% daily mos
h = 10; % horizon
for i=1:num_fold
    temp_train = k_fold_train(:,2:5,i);
    temp_test = k_fold_test(:,2:5,i);
    
    daily_train_C(:,:,:,i) = reshape_mos_daily(temp_train,h);
    daily_train_y(:,:,i) = reshape(k_fold_train(:,1,i),10,[])';
    daily_train_t(:,:,i) = reshape(k_fold_train(:,6,i),10,[])';
end

%% estimate daily beta
for i=1:num_fold
    for j = 1:h
        daily_beta(:,j,i) = daily_train_C(:,:,j,i)\daily_train_y(:,j,i);
    end
end

%% Persistent Model
% h1 = 7; % the begining hour of the day for training and test
% h2 = 16; % the last hour of the day for training and test
% tf = 13; 
% for ii=1:num_fold
%     start_test = k_fold_test_date(1,ii);
%     last_tset = k_fold_test_date(end,ii); 
% 
%     persis_test = create_persistence(start_test,last_test,measure_downsample,h1,h2);
%     persis_test.k = rearrange_persis([zeros(1,h2-h1+1); reshape(solar_test.I./solar_test.Iclr,horizon,[])'],h1,h2,tf);
%     persis_test.Ihat = persis_test.k.*reshape(solar_test.Iclr,horizon,[])';
% 
%     
% end

clrindex=[zeros(10,1,10); k_fold_test(:,1,:)./k_fold_test(:,7,:)]; % the 1st and 7th variables are measured I and Iclr
for i=1:num_fold
    clrindex_arranged(:,:,i)=rearrange_persis(reshape(clrindex(:,:,i),10,[])',7,16,13); % h1= 7:00 , h2 = 16:00, tf = 13:00
    clrindex2(:,i) = reshape(clrindex_arranged(:,:,i)',[],1);
end
persis_iclr = squeeze(k_fold_test(:,7,:)); % retrieve the 7th variable = Iclr ;
ihat_persis = clrindex2.*persis_iclr;
ihat_persis(find(ihat_persis < 0)) = 0;
ihat_persis(1:20,:) = []; % remove the initial two dates of all folds;


%% ############### lorenz #########################
%% prepare lorenz data
lorenz_train.khat=k_fold_train(:,2,:)./k_fold_train(:,7,:);
lorenz_train.C=[lorenz_train.khat.^4 k_fold_train(:,5,:).^4 lorenz_train.khat.^3 k_fold_train(:,5,:).^3 lorenz_train.khat.^2 k_fold_train(:,5,:).^2 lorenz_train.khat k_fold_train(:,5,:)];
lorenz_train.e=k_fold_train(:,2,:)-k_fold_train(:,1,:);
for i=1:num_fold
    lorenz_beta(:,i) = lorenz_train.C(:,:,i)\lorenz_train.e(:,:,i);
end
lorenz_test.khat=k_fold_test(:,2,:)./k_fold_test(:,7,:);
lorenz_test.C=[lorenz_test.khat.^4 k_fold_test(:,5,:).^4 lorenz_test.khat.^3 k_fold_test(:,5,:).^3 lorenz_test.khat.^2 k_fold_test(:,5,:).^2 lorenz_test.khat k_fold_test(:,5,:)];

%% run MOS lorenz
% predict Ihat
for i=1:num_fold
    lorenz_mos_e(:,i)=lorenz_test.C(:,:,i)*lorenz_beta(:,i);
end
ihat_lorenz=ihat_wrf-lorenz_mos_e;
ihat_lorenz(find(ihat_lorenz<0))=0;

%% ################### prepare data kf and run #######################
%% prepare data kf1a 
for fold=1:10
horizon=10;
[daily_mos_train.beta,daily_mos_train.resid,daily_mos_train.covbeta] = estimate_beta_daily(daily_train_C(:,:,:,fold),daily_train_y(:,:,fold));
% model: WRF+MOS+KF1a (measurement noise covariance is block diagonal)
h=horizon; m=h; p=size(daily_train_C,2); T=length(k_fold_train); 
daily_kf1a_test.z_init=reshape(daily_beta(:,:,fold),[],1); % beta is p x h  and z_init is ph x 1
daily_kf1a_test.A=eye(p*h);
daily_kf1a_test.W=abs(diag(daily_kf1a_test.z_init)/10000);
daily_kf1a_test.V=diag(diag(cov(daily_mos_train.resid))); % resid is D x h where D is the number of days
daily_kf1a_test.P_init=blkdiag(daily_mos_train.covbeta{1,:}); % covbeta is h-cell arrays of pxp matrices
daily_kf1a_test.C=[k_fold_test(:,2,fold) k_fold_test(:,3,fold) k_fold_test(:,4,fold) k_fold_test(:,5,fold)];
T=size(daily_kf1a_test.C,1);
daily_kf1a_test.C=reshape(mat2cell(daily_kf1a_test.C,ones(1,T), p),h,1,[]);
daily_kf1a_test.C=construct_blkdiag_kf1(daily_kf1a_test.C);
daily_kf1a_test.y=reshape(k_fold_test(:,1,fold),m,[]);
daily_kf1a_test.time=reshape(kfold_test_time,m,[]);

%% prepare data diagne 
p=3; 
diagne_test.z_init=zeros(p,1);
diagne_test.A=eye(p);
diagne_test.W=eye(p);
diagne_test.V=0.01;
diagne_test.P_init=eye(p)*5;
diagne_test.C=reshape([ones(length(k_fold_test),1) k_fold_test(:,2,fold) k_fold_test(:,5,fold)]',1,p,[]); 
diagne_test.y=(k_fold_test(:,2,fold)-k_fold_test(:,1,fold))';
diagne_test.I=k_fold_test(:,1,fold);
diagne_test.Iwrf=k_fold_test(:,2,fold);
diagne_test.time=kfold_test_time;

% adjust W to be 1e-5 x I 
diagne2_test.z_init=zeros(p,1);
diagne2_test.A=eye(p);
diagne2_test.W= 10^(-5)*eye(p);
diagne2_test.V=0.01;
diagne2_test.P_init=eye(p)*5;
diagne2_test.C=reshape([ones(length(k_fold_test),1) k_fold_test(:,2,fold) k_fold_test(:,5,fold)]',1,p,[]); 
diagne2_test.y=(k_fold_test(:,2,fold)-k_fold_test(:,1,fold))';
diagne2_test.I=k_fold_test(:,1,fold);
diagne2_test.Iwrf=k_fold_test(:,2,fold);
diagne2_test.time=kfold_test_time;


%% prepare data pelland 
p=2; 
pelland_test.z_init=zeros(p,1);
pelland_test.A=eye(p);
pelland_test.W=eye(p)*10^-5;
pelland_test.V=0.01;
pelland_test.P_init=eye(p)*5*10^-5;
pelland_test.C=reshape([ones(length(k_fold_test),1) k_fold_test(:,2,fold)]',1,p,[]); 
pelland_test.y=(k_fold_test(:,2,fold)-k_fold_test(:,1,fold))';
pelland_test.I=k_fold_test(:,1,fold);
pelland_test.Iwrf=k_fold_test(:,2,fold);
pelland_test.time=kfold_test_time;

%% ################# run kf #######################
step_hourly=30*horizon;
tf=13; % time of forecasting 
%% kf1a
daily_kf1a_test = run_kf_daily(daily_kf1a_test,tf);
%% diagne model 
diagne_test = run_kf_hourly(diagne_test,tf);
diagne_test.Ihat = reshape(diagne_test.Iwrf(horizon+1:end),horizon,[])'-diagne_test.Ihat;
diagne_test.Ihat(find(diagne_test.Ihat<0))=0;

diagne2_test = run_kf_hourly(diagne2_test,tf);
diagne2_test.Ihat = reshape(diagne2_test.Iwrf(horizon+1:end),horizon,[])'-diagne2_test.Ihat;
diagne2_test.Ihat(find(diagne2_test.Ihat<0))=0;


%% pelland model 
pelland_test = run_kf_hourly(pelland_test,tf,step_hourly);
pelland_test.Ihat = reshape(pelland_test.Iwrf(horizon+1:end),horizon,[])'-pelland_test.Ihat;
pelland_test.Ihat(find(pelland_test.Ihat<0))=0;

%% accumulate ihat of kf models
ihat_kf1a(:,fold)=reshape(daily_kf1a_test.Ihat(2:end,:)',[],1);
ihat_diagne(:,fold)=reshape(diagne_test.Ihat(2:end,:)',[],1);
ihat_diagne2(:,fold)=reshape(diagne2_test.Ihat(2:end,:)',[],1);
ihat_pelland(:,fold)=reshape(pelland_test.Ihat(2:end,:)',[],1);

end

%% ###################### ihat of all models ############
Ihat_kfold.wrf= ihat_wrf;
Ihat_kfold.persis = ihat_persis;
Ihat_kfold.mos_lorenz = ihat_lorenz;
Ihat_kfold.kf_diagne = ihat_diagne;
Ihat_kfold.kf_diagne2 = ihat_diagne2;
Ihat_kfold.kf_pelland = ihat_pelland;
Ihat_kfold.kf1a = ihat_kf1a;
Ihat_kfold.info = 'Ihat of each method (row = time, column = fold number)';

save(savefilepath,'Ihat_kfold','temp_wrf', 'p8k','measured_I') ; 

%% ################# function ######################
% function for persist
function re_persis = rearrange_persis(x,h1,h2,tf) % h1=7hrs, h2=16hrs, tf=13hrs
    p1=x(:,1:tf-h1+1);
    p2=x(:,tf-h1+2:end);
    p2=circshift(p2,1);
    p2(1,:)=0;
    re_persis = [p1 p2];
    re_persis(end,:) =[];
end


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


function specific_hour=reshape_mos_daily(C,horizon)
    for i=1:horizon
        specific_hour(:,:,i)=C(i:horizon:end,:);
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

% solve regression y=X*beta multiple times (h models)
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

%% run kf
%% hourly model
function kf_hourly=run_kf_hourly(kf_hourly,tf,varargin)
% check variable input arguments
if(length(varargin)>1)
    disp('Error: The input argument is exceed the limit'); return;
    if varargin{1}<=0 
        disp('Error: The step of adaptive noise covariance must be positive'); return;
    end
end
% predefined
A=kf_hourly.A;
C=kf_hourly.C;
W=kf_hourly.W;
V=kf_hourly.V;
zhat_time=kf_hourly.z_init; % zhat(1|0)
P_time=kf_hourly.P_init; % P(1|0)
y=kf_hourly.y;
Ihat=[];
num_sample=size(y,2);
time=kf_hourly.time;
horizon=16-7+1;
vt=[];
wt=[];
hour_left=16-tf;
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

%% daily model
function kf_daily=run_kf_daily(kf_daily,tf,varargin)
% check variable input arguments
if(length(varargin)>1)
    disp('Error: The input argument is exceed the limit'); return;
    if varargin{1}<=0
        disp('Error: The step of adaptive noise covariance must be positive'); return;
    end
end
% predefined
A=kf_daily.A;
C=kf_daily.C;
W=kf_daily.W;
V=kf_daily.V;
horizon=size(C,1);
D=size(C,3); % days = number of day
r=horizon-(max(hour(kf_daily.time(:,1)))-tf); % r = number of hours from the first hour until tf  
n=r*size(W,1)/size(V,1); 
% s=horizon-r; % s = number of hours from hour tf+1 until the end of the day
F=eye(r,horizon);
zhat_time=kf_daily.z_init;
P_time=kf_daily.P_init;
y=kf_daily.y;
time=kf_daily.time;
zhat_time=reshape(zhat_time,size(zhat_time,1),1,[]);
tf_index=find(hour(time(1:horizon))==tf); % integer index of tf
% run KF to estimate zhat(t+1|t)
for i = 1:D-1 % run from 1 to N-horizon-hour_left iterations and provide predicted value from day 2 until day D-1 
    % Measurement update at time tf
    K_tf(:,:,i) = P_time(1:n,1:n,i)*C(1:r,1:n,i)'/(C(1:r,1:n,i)*P_time(1:n,1:n,i)*C(1:r,1:n,i)'+V(1:r,1:r)); % kalman gain K(t)=P(t|t-1)C(t)'F'/F(C(t)P(t|t-1)C(t)+V(t))F' 
    zhat_mea_tf(:,i) = zhat_time(1:n,i) + K_tf(:,:,i)*(y(1:r,i)-(C(1:r,1:n,i)*zhat_time(1:n,i))); % zhat(t|t)
%     P_mea_tf(:,:,i) = P_time(:,:,i)-K(:,:,i)*F*C(:,:,i)*P_time(:,:,i); % P(t|t)
    
    % Time update at time tf
    zhat_time_tf(:,i+1) = [A(1:n,1:n)*zhat_mea_tf(:,i); zhat_time(n+1:end,i)]; % zhat(t+1|t)
%     P_time_tf(:,:,i+1) = A*P_mea_tf(:,:,i)*A' + W; % P(t+1|t)
    
    % Forecast at tf for day D+1 
    yhat(:,i+1) = C(:,:,i+1)*zhat_time_tf(:,i+1);
    
    % Measurement update 
    K(:,:,i) = P_time(:,:,i)*C(:,:,i)'/(C(:,:,i)*P_time(:,:,i)*C(:,:,i)'+V); % kalman gain K(t)=P(t|t-1)C(t)'F'/F(C(t)P(t|t-1)C(t)+V(t))F' 
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
kf_daily.K_tf=K_tf;
end