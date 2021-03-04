

%% Statistical Wilcoxon sign rank test of diffence

% H0: RMSE1 and RMSE2 have equal median

% https://www.mathworks.com/help/stats/signrank.html
% The result h is 1 if the test rejects the null hypothesis at the 5% significance level, and 0 otherwise.
% signrank(x,y,'tail','right');
% For a two-sample test, the alternate hypothesis states the data in x – y come from a distribution with median greater than 0.
% Null hypothesis : the data in x – y come from a distribution with median less than 0.
% 
% dmetric_test.rmse = zeros(num_method-1,4); % 4 columns are 1. mean of delta metric 2. h (binary) 3. p-value 4. t-stats 
% dmetrichour_test.rmse = zeros(num_method-1,4,num_hour); % 4 columns are 1. delta metric 2. h (binary) 3. p-value 4. t-stats 

% TEST_NAME = "[p,h,stats] = signrank(x,y);"; % Wilcoxin sign rank test 
TEST_NAME = "[p,h,stats] = signrank(x,y,'tail','right');"; % Wilcoxin sign rank test 

for jj = 1: (num_method -1)
    y = perfindex_kfold.rmse(:,end); % y is our method
    x = perfindex_kfold.rmse(:,jj) ; % x is each of other methods  (we know that x-y > 0 in experiments)
    eval(TEST_NAME);    
    dmetric_test.rmse(jj,:) = [mean(x-y) h p stats.signedrank];
    
    y = perfindex_kfold.mae(:,end);
    x = perfindex_kfold.mae(:,jj) ; 
    eval(TEST_NAME);
    dmetric_test.mae(jj,:) = [mean(x-y) h p stats.signedrank];
    
    y = perfindex_kfold.mbe(:,end);
    x = perfindex_kfold.mbe(:,jj) ; 
    eval(TEST_NAME);
    dmetric_test.mbe(jj,:) = [mean(x-y) h p stats.signedrank];
    
   for kk=1:num_hour
       
%        perfindex_kfoldhour.rmse is num_fold x num_method x num_hour
       y = perfindex_kfoldhour.rmse(:,end,kk); % our method
       x = perfindex_kfoldhour.rmse(:,jj,kk);
       eval(TEST_NAME);
       dmetrichour_test.rmse(jj,:,kk) = [mean(x-y) h p stats.signedrank];
       
       y = perfindex_kfoldhour.mae(:,end,kk);
       x = perfindex_kfoldhour.mae(:,jj,kk) ; 
       eval(TEST_NAME);
       dmetrichour_test.mae(jj,:,kk) = [mean(x-y) h p stats.signedrank];
       
       y = perfindex_kfoldhour.mbe(:,end,kk);
       x = perfindex_kfoldhour.mbe(:,jj,kk);
       eval(TEST_NAME);
       dmetrichour_test.mbe(jj,:,kk) = [mean(x-y) h p stats.signedrank];
   end
   
end
dmetric_test.info = {'each metric: methods x [avg(dMetric) test result (1=reject) p-value signed_rank-stats]'} ;
dmetrichour_test.info = {'each metric: methods x [avg(dMetric) test result (1=reject) p-value signed_rank-stats] x num_hour'} ;
dmetric_test
dmetrichour_test


