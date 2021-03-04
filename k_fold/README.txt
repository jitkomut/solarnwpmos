
% k-fold experiment
% Goal: 1) evaluate error with k-fold cross validation 2) use metric in each fold as each sample to run Wilcoxon test

Order of running codes

1) run_ihat.m : This produces Ihat in 10-fold
2) run_phat.m : This converts Ihat to Phat and evaluate performance
3) summary_ihat.m : Conclude performance, plot graphs, run Wilcoxon test on irradiance metrics
4) summary_phat.m : Conclude performance, plot graphs, run Wilcoxon test on solar power metrics

Required files

performance_index.m : calculate RMSE, MAE, MBE
ttest_metric.m : t-test with H0: mu_RMSE1 - mu_RMSE2 = 0
wilcoxontest_mtric.m : Wilcoxon test with H0: RMSE1-RMSE2 has a distribution with median > 0  (right tail test)
plot_metric.m : given performance indices, plot the bar graphs and print some tables
