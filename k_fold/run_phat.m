%% PV conversion
clear all;

addpath('./input_data/');
addpath('./result_files/');

load('kfold_ihat.mat')
load('pv_model.mat')

savefilepath = './result_files/kfold_phat.mat' ;

num_fold = 10; 

%% Calculate Phat
% store in Phat_kfold.XXX  as numday x num_fold array


START_TIME_INDEX = 21; % 2 days x 10 + 1 (10 = horizon)

for i=1:num_fold

Phat_kfold.wrf(:,i) = [Ihat_kfold.wrf(:,i) temp_wrf(:,i) Ihat_kfold.wrf(:,i).*temp_wrf(:,i)]*pv_conv.beta;
Phat_kfold.persis(:,i) =[Ihat_kfold.persis(:,i) temp_wrf(START_TIME_INDEX:end,i) Ihat_kfold.persis(:,i).*temp_wrf(START_TIME_INDEX:end,i)]*pv_conv.beta;
Phat_kfold.mos_lorenz(:,i) =[Ihat_kfold.mos_lorenz(:,i) temp_wrf(:,i) Ihat_kfold.mos_lorenz(:,i).*temp_wrf(:,i)]*pv_conv.beta;
Phat_kfold.kf_diagne(:,i) =[Ihat_kfold.kf_diagne(:,i) temp_wrf(START_TIME_INDEX:end,i) Ihat_kfold.kf_diagne(:,i).*temp_wrf(START_TIME_INDEX:end,i)]*pv_conv.beta;
Phat_kfold.kf_diagne2(:,i) =[Ihat_kfold.kf_diagne2(:,i) temp_wrf(START_TIME_INDEX:end,i) Ihat_kfold.kf_diagne2(:,i).*temp_wrf(START_TIME_INDEX:end,i)]*pv_conv.beta;
Phat_kfold.kf_pelland(:,i) = [Ihat_kfold.kf_pelland(:,i) temp_wrf(START_TIME_INDEX:end,i) Ihat_kfold.kf_pelland(:,i).*temp_wrf(START_TIME_INDEX:end,i)]*pv_conv.beta;
Phat_kfold.kf1a(:,i) = [Ihat_kfold.kf1a(:,i) temp_wrf(START_TIME_INDEX:end,i) Ihat_kfold.kf1a(:,i).*temp_wrf(START_TIME_INDEX:end,i)]*pv_conv.beta;

%% Threshold negative power to zero

Phat_kfold.wrf(find( Phat_kfold.wrf <0)) = 0;
Phat_kfold.persis(find(Phat_kfold.persis < 0))=0;
Phat_kfold.mos_lorenz(find(Phat_kfold.mos_lorenz < 0 ))=0;
Phat_kfold.kf_diagne(find( Phat_kfold.kf_diagne <0))=0;
Phat_kfold.kf_diagne2(find( Phat_kfold.kf_diagne2 <0))=0;
Phat_kfold.kf_pelland (find(Phat_kfold.kf_pelland < 0 ))=0;
Phat_kfold.kf1a  (find(Phat_kfold.kf1a < 0 ))=0;

end

%% Performance Evaluation
% START_TIME_INDEX = 11; % Supachai original value (Lagged time for KF
% models)
INSTALL_CAP = 8;

num_method = 7 ; num_folds = 10; num_hour = 10;
method_label = {'wrf','persis','mos_lorenz','kf_pelland','kf_diagne','kf_diagne2','kf1a'};
method_printlabel = {'WRF','Persistence','Mos_{Lorenz}','KF_{Pelland}','KF_{Diagne}','KF_{Diagne2}','KF_{daily}'};

% Metric in each fold
perfindex_kfold.rmse = zeros(num_folds,num_method);
perfindex_kfold.mae = zeros(num_folds,num_method);
perfindex_kfold.mbe = zeros(num_folds,num_method);

% Metric in each fold and each hour
perfindex_kfoldhour.rmse = zeros(num_folds,num_method,num_hour);
perfindex_kfoldhour.mae = zeros(num_folds,num_method,num_hour);
perfindex_kfoldhour.mbe = zeros(num_folds,num_method,num_hour);

for ii=1:num_folds
    for jj=1:num_method
    
        % Evaluate the performance index
        % for KF methdos, evaluation start at index 21
        
        if strcmp(method_label{jj},'wrf') || strcmp(method_label{jj},'mos_lorenz') 
            var_name = ['performance_index(Phat_kfold.',method_label{jj},'(:,ii),p8k(:,ii))'];
        else
            var_name = ['performance_index(Phat_kfold.',method_label{jj},'(:,ii),p8k(START_TIME_INDEX:end,ii))'];
        end
        
        [perf,perfhour] = eval(var_name);
        
        perfindex_kfold.rmse(ii,jj) = perf.rmse*100/INSTALL_CAP;
        perfindex_kfoldhour.rmse(ii,jj,:) = perfhour.rmse'*100/INSTALL_CAP; 
        
        perfindex_kfold.mae(ii,jj) = perf.mae*100/INSTALL_CAP;
        perfindex_kfoldhour.mae(ii,jj,:) = perfhour.mae'*100/INSTALL_CAP; 
        
        perfindex_kfold.mbe(ii,jj) = perf.mbe*100/INSTALL_CAP;
        perfindex_kfoldhour.mbe(ii,jj,:) = perfhour.mbe'*100/INSTALL_CAP; 
    end
end
perfindex_kfold.info = {'Unit of solar power is normalized by capacity (%)'};
perfindex_kfoldhour.info = {'Unit of solar power is normalized by capacity (%)'};

save(savefilepath,'perfindex_kfold', 'perfindex_kfoldhour','method_label', 'method_printlabel');
