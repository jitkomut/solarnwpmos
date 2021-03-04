
%% Plot graphs

% method_printlabel = {'WRF','Mos_{Lorenz}','KF_{Diagne}','KF_{Pelland}','KF_{daily}'};
hour_label = {'7:00','8:00','9:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00'};

%% overall metric

figure(1); fig = tiledlayout(3,1); fig.TileSpacing = 'compact'; fig.Padding = 'compact';

nexttile;
bar(mean(perfindex_kfold.rmse,1)); % average over k-folds
ylabel('RMSE');
set(gca,'XTickLabel',method_printlabel,'fontsize',28);


nexttile;
bar(mean(perfindex_kfold.mae,1));
ylabel('MAE');
set(gca,'XTickLabel',method_printlabel,'fontsize',28);

nexttile;
bar(mean(perfindex_kfold.mbe,1));
ylabel('MBE');
set(gca,'XTickLabel',method_printlabel,'fontsize',28);

%% Print overall metric 

M = [mean(perfindex_kfold.rmse,1) ; mean(perfindex_kfold.mae,1) ; mean(perfindex_kfold.mbe,1)]'; % methods x [RMSE MAE MBE]

row_head = {'WRF','Persis','Lorenz','Pelland','Diagne','Diagne_adj','KF_daily'};
col_head = {'Methods','RMSE','MAE','MBE'};

disp('Overall metric table');
printtable_head(M,row_head,col_head,'%.2f');


%% Performance of each hour

figure(2); fig = tiledlayout(3,1); fig.TileSpacing = 'compact'; fig.Padding = 'compact';

% avg_metrichour = method x hour
avg_metrichour.rmse = squeeze(mean(perfindex_kfoldhour.rmse,1))'; % W/m^2
avg_metrichour.mae = squeeze(mean(perfindex_kfoldhour.mae,1))'; 
avg_metrichour.mbe = squeeze(mean(perfindex_kfoldhour.mbe,1))'; 

nexttile;

bar(avg_metrichour.rmse); % average over k-folds
grid on; ylabel('RMSE');
legend(method_printlabel,'location','eastoutside','fontsize',20);
set(gca,'XTickLabel',hour_label,'fontsize',28);

nexttile;
bar(avg_metrichour.mae);
grid on; ylabel('MAE');
set(gca,'XTickLabel',hour_label,'fontsize',28);

nexttile;
bar(avg_metrichour.mbe);
grid on; ylabel('MBE');
set(gca,'XTickLabel',hour_label,'fontsize',28);
