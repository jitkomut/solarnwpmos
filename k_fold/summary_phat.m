
clear all ;

addpath('./result_files/');

load('kfold_phat.mat');

num_method = length(method_label); 
num_folds = 10; num_hour = 10;

% method_label = {'wrf','mos_lorenz','kf_diagne','kf_pelland','kf1a'};


%% Plot graphs

% method_printlabel = {'WRF','Mos_{Lorenz}','KF_{Diagne}','KF_{Pelland}','KF_{daily}'};

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

disp('Print statistic and p-value');
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

