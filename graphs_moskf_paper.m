
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

% This file plots all the graphs and type tables in the paper.

%% Irradiance performance
clear all; close all;
addpath('./result_files/');

load moskf_test_results
% load moskf_train_results

% Test dates are Jul 2018 - Dec 2018

method_name = {'wrf','persis','lorenz','diagne','pelland','daily_kf1a'};

label = {'WRF','Persistence','Mos_{Lorenz}','KF_{Diagne}','KF_{Pelland}','KF_{daily}'};

hour_label = {'7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00'};


% matrix of performance index
rmse_i = []; mbe_i = []; rmse_i_hour = [];  mbe_i_hour = [];
for k = 1:length(method_name)
    name = [method_name{k} '_test.perf.rmse'];
    rmse_i = [rmse_i eval(name)];
    
    name = [method_name{k} '_test.perf.mbe'];
    mbe_i = [mbe_i eval(name)];
    
    name = [method_name{k} '_test.perf_specific.rmse'];
    rmse_i_hour = [rmse_i_hour ; eval(name)];
    
    name = [method_name{k} '_test.perf_specific.mbe'];
    mbe_i_hour = [mbe_i_hour ; eval(name)];
end


% Performance of Irradiance forecasting (by each hour)
t = tiledlayout(2,1,'TileSpacing','Compact');

nexttile;
bar(1000*rmse_i_hour'); ylabel('RMSE of irradiance (W/sqm)','fontsize',20);
legend(label,'location','northwest','fontsize',18);
set(gca,'XTickLabel',hour_label,'fontsize',20);

nexttile;
bar(1000*mbe_i_hour'); ylabel('MBE of irradiance (W/sqm)','fontsize',20);
xlabel('Hours','fontsize',17);
set(gca,'XTickLabel',hour_label,'fontsize',20);
set(gcf,'WindowState','fullscreen')

figfilename = 'perf_irradince';
% print(figfilename,'-painters','-depsc','-r300');

%% Solar power performance
clear all; close all

SOLAR_PLANT = 1 ; % 8kW
% SOLAR_PLANT = 2 ; % 15kW

if SOLAR_PLANT == 1
    load moskf_conversion8kW_test_results
    install_cap = 8;
    figfilename = 'perf_power8kW';
elseif SOLAR_PLANT == 2
    load moskf_conversion15kW_test_results
    install_cap = 15;
    figfilename = 'perf_power15kW';
end

method_name = {'wrf','persis','lorenz','diagne','pelland','daily_kf1a'};
label = {'WRF','Persistence','Mos_{Lorenz}','KF_{Diagne}','KF_{Pelland}','KF_{daily}'};
hour_label = {'7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00'};

% matrix of performance index
rmse_p = []; mbe_p = []; % num_hour x 1
rmse_p_hour = []; nrmse_p_hour = [];  mbe_p_hour = []; nmbe_p_hour = [];% num_method x num_hour
for k = 1:length(method_name)
    name = [method_name{k} '_test.perf_p.rmse'];
    rmse_p = [rmse_p eval(name)];
    
    name = [method_name{k} '_test.perf_p.mbe'];
    mbe_p = [mbe_p eval(name)];
    
    name = [method_name{k} '_test.perf_p_specific.rmse'];
    rmse_p_hour = [rmse_p_hour ; eval(name)];
    
    name = [method_name{k} '_test.perf_p_specific.nrmse_cap'];
    nrmse_p_hour = [nrmse_p_hour ; eval(name)];
    
    name = [method_name{k} '_test.perf_p_specific.mbe'];
    mbe_p_hour = [mbe_p_hour ; eval(name)];
    
    name = [method_name{k} '_test.perf_p_specific.nmbe_cap'];
    nmbe_p_hour = [nmbe_p_hour ; eval(name)];
    
end


% Performance of Irradiance forecasting (by each hour)
t = tiledlayout(2,1,'TileSpacing','Compact');

nexttile;
bar(nrmse_p_hour'); ylabel('NRMSE of solar power (%)','fontsize',20);
title([num2str(install_cap), '-kW Solar Station']);
legend(label,'location','northwest','fontsize',16);
set(gca,'XTickLabel',hour_label,'fontsize',20)

nexttile;
bar(nmbe_p_hour'); ylabel('NMBE of solar power (%)','fontsize',20);
xlabel('Hours','fontsize',20);
set(gca,'XTickLabel',hour_label,'fontsize',20);
set(gcf,'WindowState','fullscreen')


% print(figfilename,'-painters','-depsc','-r300');

%% Time series plot
clear all; close all; clc;

[method_label,rankI] = rankerr ; 

t = tiledlayout(2,1,'TileSpacing','Compact');
nexttile;
plot(1000*rankI.bestmea,'*-','markersize',8); hold on;
plot(1000*rankI.bestval,'.--','markersize',14); 
legend(method_label,'location','northwest','fontsize',16);
title('10-day lowest daily-averaged RMSE','fontsize',20);
ylabel('Irradiance (W/sqm)','fontsize',20); %xlabel('Date','fontsize',15);
set(gca,'XTickLabel',rankI.besttime,'fontsize',20,'YLim',[-10 1100]);

nexttile;
plot(1000*rankI.worstmea,'*-','markersize',8); hold on;
plot(1000*rankI.worstval,'.--','markersize',14); %legend(method_label,'location','northwest');
ylabel('Irradiance (W/sqm)','fontsize',20); %xlabel('Date','fontsize',15);
title('10-day highest daily-averaged RMSE','fontsize',20)
set(gca,'XTickLabel',rankI.worsttime,'fontsize',20,'YLim',[-10 1100]);
set(gcf,'WindowState','fullscreen')

figfilename = 'timeseries_I'
% print(figfilename,'-painters','-depsc','-r300');

%% Result Tables

clear all; close all;

load moskf_test_results
% load moskf_train_results

% Test dates are Jul 2018 - Dec 2018

method_name = {'wrf','persis','lorenz','diagne','pelland','daily_kf1a'};
% Our method name is 'daily_kf1a'

% matrix of performance index
rmse_i = []; mbe_i = []; 

for k = 1:length(method_name)
    name = [method_name{k} '_test.perf.rmse'];
    rmse_i = [rmse_i eval(name)];
    
    name = [method_name{k} '_test.perf.mbe'];
    mbe_i = [mbe_i eval(name)];
        
end
avg_i = mean(daily_kf1a_test.y(:)); % averaged GHI
nrmse_i = 100*rmse_i/avg_i; nmbe_i = 100*mbe_i/avg_i  ; % percent

clearvars -except rmse_i mbe_i nrmse_i nmbe_i method_name

load moskf_conversion8kW_test_results

nrmse_p8kW = []; nmbe_p8kW = []; 

for k = 1:length(method_name)
    name = [method_name{k} '_test.perf_p.nrmse_cap'];
    nrmse_p8kW = [nrmse_p8kW eval(name)];
    
    name = [method_name{k} '_test.perf_p.nmbe_cap'];
    nmbe_p8kW = [nmbe_p8kW eval(name)];
        
end

clearvars -except rmse_i mbe_i nrmse_i nmbe_i  nrmse_p8kW nmbe_p8kW method_name

load moskf_conversion15kW_test_results

nrmse_p15kW = []; nmbe_p15kW = []; 

for k = 1:length(method_name)
    name = [method_name{k} '_test.perf_p.nrmse_cap'];
    nrmse_p15kW = [nrmse_p15kW eval(name)];
    
    name = [method_name{k} '_test.perf_p.nmbe_cap'];
    nmbe_p15kW = [nmbe_p15kW eval(name)];
        
end

clearvars -except rmse_i mbe_i nrmse_i nmbe_i  nrmse_p8kW nmbe_p8kW nrmse_p15kW nmbe_p15kW method_name

M = [1000*rmse_i ; 1000*nrmse_i ; nrmse_p8kW ; nrmse_p15kW ; 1000*mbe_i ; nmbe_p8kW ; nmbe_p15kW];
[m,n] = size(M);
method_label = {'WRF','Persistent','Lorenz','Diagne','Pelland','KF (daily)'};
for k=1:n
    bold_table_head{k} = [' \\bf ' method_label{k}]; % all latex command starting with \ need to be typed twice
end
printtable(M,[],bold_table_head)

M1 = [1000*rmse_i ; nrmse_i ;  nrmse_p8kW ; nrmse_p15kW ]'; col_head1  = {'Methods','RMSE (W/sqm)','NRMSE','NRMSE (8kW)','NRMSE (15kW)'};
M2 = [1000*mbe_i ; nmbe_i ; nmbe_p8kW ; nmbe_p15kW]'; col_head2 = {'Methods','MBE (W/sqm)','NMBE','NMBE (8kW)','NMBE (15kW)'};

printtable_head(M1,method_label,col_head1);
printtable_head(M2,method_label,col_head2);




