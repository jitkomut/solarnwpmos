
function[method_label,rankI] = rankerror


% Validation set 1
% Test set: Jul 2018 - Dec 2018
% Pelland, Diagne, KF Ihat starts from Jul 2 in test set. 
% Measurement and Lorenz starts from Jul 1 in test set. 
% calculate error on Jul 2, 7:00 to Dec 31, 16:00 (183 days)
% (originally, saved in daresult_validation1.mat)

% Validation set 2
% Test set: Jan 2018 - Dec 2018
% Pelland, Diagne, KF Ihat starts from Jan 2 in test set. 
% Measurement and Lorenz starts from Jan 1 in test set. 
% calculate error on Jan 2, 7:00 to Dec 31, 16:00 (364 days)
% (originally, saved in daresult_validation2.mat)

addpath('./result_files/');

load moskf_test_results
load solar_splitdata

num_hour = 10; num_days = 183 ; % hard code
numday2plot = 10;

method_name = {'lorenz_test','diagne_test','pelland_test','daily_kf1a_test'};
method_label = {'Measurement','MOS (lorenz)','KF (diagne)','KF (pelland)','KF (daily)'};

datetimevec = reshape(solar_test.time(num_hour+1:end),num_hour,num_days)'; % in the format of num_days x num_hour 
% datetimevec = daily_kf1a_test.time ; % num_days x num_hour. THis agrees
% with line 128

%% Calculate the best and worst predicted Irradiance
Imea = solar_test.I(num_hour+1:end); 
Imea = reshape(Imea,num_hour,num_days)'; % num_days x num_hour

err = zeros(num_days,num_hour,length(method_name)); % num_days x num_hour x methods
for k=1:length(method_name)
    if strcmp(method_name{k},'lorenz_test')
        err(:,:,k) = eval([method_name{k} '.Ihat(2:end,:)'])-Imea ; % start from Jul 2/Jan 2 (index 2)
    else
        err(:,:,k) = eval([method_name{k} '.Ihat'])-Imea ; % start from Jul 2/Jan 2 (index 1)
    end
end

RMSE_allmethod = squeeze(sqrt(sum(sum(err.^2,1),2)/(num_days*num_hour)))

avgRMSE1day = squeeze(sqrt(sum(err.^2,2)/num_hour)); % 2-norm of error for each day and each method, size = num_days x methods
[rankRMSE,rankday_ind] = sort(avgRMSE1day(:,end),'ascend'); % rank the best days of our method
% [rankRMSE,rankday_ind] = sort(avgRMSE1day(:,1),'ascend'); % rank the best days of Lorenz
% [rankRMSE,rankday_ind] = sort(avgRMSE1day(:,2),'ascend'); % rank the best days of Diagne
% [rankRMSE,rankday_ind] = sort(avgRMSE1day(:,3),'ascend'); % rank the best days of Pelland

bestday_ind = datetimevec(rankday_ind(1:numday2plot),:) ; 
bestday2plot = zeros(numday2plot*num_hour,length(method_name)); % to store time series of best days

worstday_ind = datetimevec(rankday_ind(end:-1:end-numday2plot+1),:) ; 
worstday2plot = zeros(numday2plot*num_hour,length(method_name)); % to store time series of best days

for k=1:length(method_name)
    [C,IA,IB] = intersect(bestday_ind(:,1),eval([method_name{k} '.time(:,1)']),'rows','stable'); % compare only dates at 7:00 (index 1)
    tmp = eval([method_name{k} '.Ihat(IB,:)']); % select the time series of best day , numday2plot x numhour 
    bestday2plot(:,k) = reshape(tmp',numday2plot*num_hour,1);
    
    [C,IA,IB] = intersect(worstday_ind(:,1),eval([method_name{k} '.time(:,1)']),'rows','stable'); % compare only dates at 7:00 (index 1)
    tmp = eval([method_name{k} '.Ihat(IB,:)']); % select the time series of worst day , numday2plot x numhour 
    worstday2plot(:,k) = reshape(tmp',numday2plot*num_hour,1);
end
[C,IA,IB] = intersect(bestday_ind(:,1),datetimevec(:,1),'rows','stable'); % compare only dates at 7:00 (index 1)
tmp = Imea(IB,:); bestImea = reshape(tmp',numday2plot*num_hour,1);

[C,IA,IB] = intersect(worstday_ind(:,1),datetimevec(:,1),'rows','stable'); % compare only dates at 7:00 (index 1)
tmp = Imea(IB,:); worstImea = reshape(tmp',numday2plot*num_hour,1);

bestday_label = datestr(bestday_ind(:,1),'dd-mmm-yy'); 
worstday_label = datestr(worstday_ind(:,1),'dd-mmm-yy'); 

rankI.bestmea = bestImea ;
rankI.bestval = [bestday2plot];
rankI.besttime = bestday_label;

rankI.worstmea = worstImea ;
rankI.worstval = [worstday2plot];
rankI.worsttime = worstday_label;

