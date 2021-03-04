
clear all
% load('ihat_kfold.mat')

addpath('./result_files/');

load('kfold_ihat.mat');

method_label = {'wrf','persis','mos_lorenz','kf_pelland','kf_diagne','kf_diagne2','kf1a'};

num_method = length(method_label) ; num_folds = 10; num_hour = 10;

% Metric in each fold
perfindex_kfold.rmse = zeros(num_folds,num_method);
perfindex_kfold.mae = zeros(num_folds,num_method);
perfindex_kfold.mbe = zeros(num_folds,num_method);

% Metric in each fold and each hour
perfindex_kfoldhour.rmse = zeros(num_folds,num_method,num_hour);
perfindex_kfoldhour.mae = zeros(num_folds,num_method,num_hour);
perfindex_kfoldhour.mbe = zeros(num_folds,num_method,num_hour);

%% calculate performance indices
for ii=1:num_folds
    for jj=1:num_method
    
        % Evaluate the performance index
        % for KF methdos and persistence, evaluation start at index 21 (2days x 10 + 1)
        if strcmp(method_label{jj},'wrf') || strcmp(method_label{jj},'mos_lorenz') 
            var_name = ['performance_index(Ihat_kfold.',method_label{jj},'(:,ii),measured_I(:,ii))'];
        else
            var_name = ['performance_index(Ihat_kfold.',method_label{jj},'(:,ii),measured_I(21:end,ii))'];
        end
        
        [perf,perfhour] = eval(var_name);
        
        perfindex_kfold.rmse(ii,jj) = 1000*perf.rmse;
        perfindex_kfoldhour.rmse(ii,jj,:) = 1000*perfhour.rmse'; 
        
        perfindex_kfold.mae(ii,jj) = 1000*perf.mae;
        perfindex_kfoldhour.mae(ii,jj,:) = 1000*perfhour.mae'; 
        
        perfindex_kfold.mbe(ii,jj) = perf.mbe;
        perfindex_kfoldhour.mbe(ii,jj,:) = 1000*perfhour.mbe'; 
    end
end
perfindex_kfold.info = {'Unit of irradiance is W/sqm'};
perfindex_kfoldhour.info = {'Unit of irradiance is W/sqm'};

%% Plot graphs

method_printlabel = {'WRF','Persistence','Mos_{Lorenz}','KF_{Pelland}','KF_{Diagne}','KF_{Diagne2}','KF_{daily}'};

plot_metric;


%% t-test

ttest_metric ;

disp('t-test of RMSE: dRMSE & test result (1= reject null) & p-value & t-stat')
dmetric_test.rmse

disp('t-test of MAE: dMAE & test result (1= reject null) & p-value & t-stat')
dmetric_test.mae

disp('t-test of MBE: dMBE & test result (1= reject null) & p-value & t-stat')
dmetric_test.mbe

%% Wilcoxon signed rank test

wilcoxontest_metric ;

disp('Wilcoxon test of RMSE: dRMSE & test result (1= reject null) & p-value & signed_rank')
dmetric_test.rmse

disp('Wilcoxon test of MAE: dMAE & test result (1= reject null) & p-value & signed_rank')
dmetric_test.mae

disp('Wilcoxon test of MBE: dMBE & test result (1= reject null) & p-value & signed_rank')
dmetric_test.mbe


M = [dmetric_test.rmse(:,3) dmetric_test.mae(:,3) dmetric_test.mbe(:,3)]; 
BFcode = [dmetric_test.rmse(:,2) dmetric_test.mae(:,2) dmetric_test.mbe(:,2)]; 
BFcode = (BFcode == 1); % convert to logical

row_head = {'WRF','Persis','Lorenz','Pelland','Diagne','Diagne (adjusted)','KF_daily'};
col_head = {'Metrics','RMSE','MAE','MBE'};

disp('Print only p-values');
printtable_head(M,row_head,col_head,'%.2e',BFcode);



M = [dmetric_test.rmse(:,[4 3]) dmetric_test.mae(:,[4 3]) dmetric_test.mbe(:,[4 3])]; 
BFcode = [zeros(num_method-1,1) dmetric_test.rmse(:,2) zeros(num_method-1,1) dmetric_test.mae(:,2) zeros(num_method-1,1) dmetric_test.mbe(:,2)]; 
BFcode = (BFcode == 1); % convert to logical

printFORMAT = cell(1,3*2); % 3 metrics x 2
for jj=1: 3*2 
    if mod(jj,2)
        printFORMAT{jj} = '%d';
    else
        printFORMAT{jj} = '%.2e';
    end
end

col_head = {'Metrics'}; 
for jj=1:3  % 3 metrics
    col_head = [col_head  {'stat','p-val'}];
end


printtable_head(M,row_head,col_head,printFORMAT,BFcode);

