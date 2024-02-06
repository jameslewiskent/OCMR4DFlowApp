% function     stats = StatsFlow(stasis_atria,stasis_atria_erode,stasis_atria_wall,peak_atria,flow_mag_atria,vox_vol,time_stamps,vorticity_mag_atria,vortex_sz_ml,lambda2_atria,vorticity_map_mask,vorticity_mag_wall,vorticity_mag_internal,peak_internal,peak_wall,flow_mag_internal,flow_mag_wall,vorticity_map_internal, vorticity_map_wall);
function stats = StatsFlowOriginal(stasis_atria,peak_atria,flow_mag_atria,vox_vol,time_stamps,vorticity_mag_atria,vortex_sz_ml,lambda2_atria);
%%for time analysis add ,stasis_atria_systole,stasis_atria_diastole,stasis_atria_systole_erode,stasis_atria_systole_wall,stasis_atria_diastole_erode,stasis_atria_diastole_wall,peak_atria_systole,peak_internal_systole,peak_wall_systole,peak_atria_diastole,peak_internal_diastole,peak_wall_diastole

     stats.stasis_min_max = [min(stasis_atria),max(stasis_atria)];
    stats.stasis_iqr = iqr(stasis_atria);
    stats.stasis_median = median(stasis_atria);
    stats.stasis_mean = mean(stasis_atria);
    stats.stasis_std = std(stasis_atria);
    stats.stasis_vol_0 = sum(stasis_atria>0)*vox_vol;
    stats.stasis_vol_10 = sum(stasis_atria>0.1)*vox_vol;
    stats.stasis_vol_20 = sum(stasis_atria>0.2)*vox_vol;
    stats.stasis_vol_30 = sum(stasis_atria>0.3)*vox_vol;
    stats.stasis_vol_40 = sum(stasis_atria>0.4)*vox_vol;
    stats.stasis_vol_50 = sum(stasis_atria>0.5)*vox_vol;
    stats.stasis_vol_60 = sum(stasis_atria>0.6)*vox_vol;
    stats.stasis_vol_70 = sum(stasis_atria>0.7)*vox_vol;
    stats.stasis_vol_80 = sum(stasis_atria>0.8)*vox_vol;
    stats.stasis_vol_90 = sum(stasis_atria>0.9)*vox_vol;
    stats.stasis_vol_100 = sum(stasis_atria>1)*vox_vol;

    stats.peak_min_max = [min(peak_atria),max(peak_atria)];
    stats.peak_iqr = iqr(peak_atria);
    stats.peak_median = median(peak_atria);
    stats.peak_mean = mean(peak_atria);

    stats.peak_std = std(peak_atria);

    stats.mean_flow = mean(flow_mag_atria(:));

    stats.mean_flow_auc = sum(mean(flow_mag_atria,1))*mean(diff(time_stamps));

    stats.quality_factor_std_vel = sqrt(mean(var(diff(flow_mag_atria,1,2),[],2),1));
    
%     stats.vorticity_time = mean(vorticity_atria,1);
%     stats.cov = mean(cov_atria,1);
    stats.vorticity_int = sum(mean(vorticity_mag_atria,1))*mean(diff(time_stamps));
%   
    stats.vorticity_peak = max(mean(vorticity_mag_atria,1));

    stats.vortexSize_peak = max(vortex_sz_ml);
    stats.vortexsize_int = sum((mean(vortex_sz_ml,1))*mean(diff(time_stamps)));
    stats.vortexsize_mean = mean(vortex_sz_ml);
    stats.vortexSize_change = max(vortex_sz_ml)-min(vortex_sz_ml);
    stats.lambda2_mean = mean(lambda2_atria(:));
    stats.lambda2_min_mean = min(mean(lambda2_atria));
    stats.lambda2_threshold = mean(lambda2_atria(:))/2;
%     stats.cov_int = sum(mean(cov_atria,1));
%     stats.intersection_lowvort_histas_wall = (sum(vorticity_map_wall<15&stasis_atria_wall>0.7)*vox_vol)/(sum(stasis_atria_wall>0.7)*vox_vol)*100;
end