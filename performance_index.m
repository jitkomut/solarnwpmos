%% Performance indices
function [performance,performance_hour]=performance_index(predicted,measure)
predicted=reshape(predicted',[],1);
% overall
performance.rmse=sqrt(immse(predicted,measure)); 
performance.mae=mean(abs(predicted-measure)); 
performance.mape=mean(abs((predicted-measure)./measure))*100;
performance.mbe=mean(predicted-measure); 
performance.nrmse_mean=performance.rmse/mean(measure)*100;
performance.nrmse_sd=performance.rmse/std(measure)*100;
performance.nrmse_maxmin=performance.rmse/(max(measure)-min(measure))*100;
% specified hour
error=predicted-measure;
error=reshape(error,10,[])';
for i=1:1:10
    performance_hour.rmse(i)=norm(error(:,i))/sqrt(length(error(:,i)));
    performance_hour.mae(i)=mean(abs(error(:,i))); 
    performance_hour.mape(i)=mean(abs(error(:,i)./measure(i:10:end)))*100;
    performance_hour.mbe(i)=mean(error(:,i)); 
    performance_hour.nrmse_mean(i)=performance.rmse/mean(measure(i:10:end))*100;
    performance_hour.nrmse_sd(i)=performance.rmse/std(measure(i:10:end))*100;
    performance_hour.nrmse_maxmin(i)=performance.rmse/(max(measure(i:10:end))-min(measure(i:10:end)))*100;
end
end
