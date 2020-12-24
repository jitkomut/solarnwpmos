%% Performance indices
function [performance,performance_hour]=performance_index_pv(predicted,measure,install_cap)
% overall
performance.rmse=sqrt(immse(predicted,measure)); 
performance.mae=mean(abs(predicted-measure)); 

    tmp = abs((predicted-measure)./measure); 
    tmp(isnan(tmp)) = []; % remove NaN from division by zero
    tmp(isinf(tmp)) = []; % remove Inf from division by zero
  
performance.mape=mean(tmp)*100;


performance.mbe=mean(predicted-measure); 
performance.nrmse_cap = performance.rmse*100/install_cap;
performance.nrmse_mean=performance.rmse*100/mean(measure);
performance.nrmse_sd=performance.rmse*100/std(measure);
performance.nrmse_maxmin=performance.rmse*100/(max(measure)-min(measure));
performance.nmbe_cap = performance.mbe*100/install_cap;

% specified hour
error=predicted-measure;
error=reshape(error,10,[])';
for i=1:1:10
    performance_hour.rmse(i)=norm(error(:,i))/sqrt(length(error(:,i)));
    performance_hour.mae(i)=mean(abs(error(:,i))); 
    
    tmp = abs(error(:,i)./measure(i:10:end)); 
    tmp(isnan(tmp)) = []; % remove NaN from division by zero
    tmp(isinf(tmp)) = []; % remove Inf from division by zero
    performance_hour.mape(i)=mean(tmp)*100;
    performance_hour.mbe(i)=mean(error(:,i)); 
end

performance_hour.nrmse_cap = performance_hour.rmse*100/install_cap;
performance_hour.nrmse_mean = performance_hour.rmse*100/mean(measure(:));
performance_hour.nrmse_sd = performance_hour.rmse*100/std(measure(:));
performance_hour.nrmse_maxmin = performance_hour.rmse*100/(max(measure(:))-min(measure(:)));

performance_hour.nmbe_cap = performance_hour.mbe*100/install_cap;
end
