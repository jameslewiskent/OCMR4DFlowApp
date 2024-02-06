function     stats = StatsFlowThesis(stasis_atria,stasis_atria_erode,stasis_atria_wall,peak_atria,flow_mag_atria,vox_vol,time_stamps,vorticity_mag_atria,vortex_sz_ml,lambda2_atria,vorticity_map_mask,vorticity_mag_wall,vorticity_mag_internal,peak_internal,peak_wall,flow_mag_internal,flow_mag_wall,vorticity_map_internal, vorticity_map_wall)

%%for time analysis add ,stasis_atria_systole,stasis_atria_diastole,stasis_atria_systole_erode,stasis_atria_systole_wall,stasis_atria_diastole_erode,stasis_atria_diastole_wall,peak_atria_systole,peak_internal_systole,peak_wall_systole,peak_atria_diastole,peak_internal_diastole,peak_wall_diastole

    
    stats.stasis_systole = mean(stasis_atria_systole);
    stats.stasis_diastole = mean(stasis_atria_diastole);
   
    stats.peak_systole = mean(peak_atria_systole);
    stats.peak_diastole = mean(peak_atria_diastole);

    stats.mean_flow_systole = mean(flow_mag_atria_systole(:));
    stats.mean_flow_diastole = mean(flow_mag_atria_diastole(:));

    stats.vorticity_int_systole = sum(mean(vorticity_mag_atria_systole,1))*mean(diff(time_stamps));
    stats.vorticity_int_diastole = sum(mean(vorticity_mag_atria_diastole,1))*mean(diff(time_stamps));
    stats.vorticity_peak_systole = max(mean(vorticity_mag_atria_systole,1));
    stats.vorticity_peak_diastole = max(mean(vorticity_mag_atria_diastole,1));

    

%%%%%%%%%%%%%%%%%%% INTERNAL %%%%%%%%%%%%%%%%%%%%%%
    stats.stasis_internal_min_max = [min(stasis_atria_erode),max(stasis_atria_erode)];
    stats.stasis_internal_iqr = iqr(stasis_atria_erode);
    stats.stasis_internal_median = median(stasis_atria_erode);
    stats.stasis_internal_mean = mean(stasis_atria_erode);
    stats.stasis_systole_internal=mean(stasis_atria_systole_erode);
    stats.stasis_diastole_internal=mean(stasis_atria_diastole_erode);
    stats.stasis_internal_std = std(stasis_atria_erode);

    stats.peak_internal_min_max = [min(peak_internal),max(peak_internal)];
    stats.peak_internal_iqr = iqr(peak_internal);
    stats.peak_internal_median = median(peak_internal);
    stats.peak_internal_mean = mean(peak_internal);
    stats.peak_internal_systole_mean = mean(peak_internal_systole);
    stats.peak_internal_diastole_mean = mean(peak_internal_diastole);
    stats.peak_internal_std = std(peak_internal);

    stats.mean_internal_flow = mean(flow_mag_internal(:));
    stats.mean_internal_flow_systole = mean(flow_mag_internal_systole(:));
    stats.mean_internal_flow_diastole = mean(flow_mag_internal_diastole(:));

    stats.vorticity_internal_int = sum(mean(vorticity_mag_internal,1))*mean(diff(time_stamps));
    stats.vorticity_internal_int_systole = sum(mean(vorticity_mag_internal_systole,1))*mean(diff(time_stamps));
    stats.vorticity_internal_int_diastole = sum(mean(vorticity_mag_internal_diastole,1))*mean(diff(time_stamps));

    stats.vorticity_internal_peak = max(mean(vorticity_mag_internal,1));
    stats.vorticity_internal_peak_systole = max(mean(vorticity_mag_internal_systole,1));
    stats.vorticity_internal_peak_diastole = max(mean(vorticity_mag_internal_diastole,1));


%%%%%%%%%%%%%%% WALL %%%%%%%%%%%%%%

    stats.stasis_wall_min_max = [min(stasis_atria_wall),max(stasis_atria_wall)];
    stats.stasis_wall_iqr = iqr(stasis_atria_wall);
    stats.stasis_wall_median = median(stasis_atria_wall);
    stats.stasis_wall_mean = mean(stasis_atria_wall);
    stats.stasis_systole_wall=mean(stasis_atria_systole_wall);
    stats.stasis_diastole_wall=mean(stasis_atria_diastole_wall);
    stats.stasis_wall_std = std(stasis_atria_wall);

    stats.peak_wall_min_max = [min(peak_wall),max(peak_wall)];
    stats.peak_wall_iqr = iqr(peak_wall);
    stats.peak_wall_median = median(peak_wall);
    stats.peak_wall_mean = mean(peak_wall);
    stats.peak_wall_systole_mean = mean(peak_internal_systole);
    stats.peak_wall_diastole_mean = mean(peak_internal_diastole);
    stats.peak_wall_std = std(peak_wall);

    stats.mean_wall_flow = mean(flow_mag_wall(:));
    stats.mean_wall_flow_systole = mean(flow_mag_wall_systole(:));
    stats.mean_wall_flow_diastole = mean(flow_mag_wall_diastole(:));

    stats.vorticity_wall_int = sum(mean(vorticity_mag_wall,1))*mean(diff(time_stamps));
    stats.vorticity_wall_int_systole = sum(mean(vorticity_mag_wall_systole,1))*mean(diff(time_stamps));
    stats.vorticity_wall_int_diastole = sum(mean(vorticity_mag_wall_diastole,1))*mean(diff(time_stamps));
    stats.vorticity_wall_peak = max(mean(vorticity_mag_wall,1));
    stats.vorticity_wall_peak_systole = max(mean(vorticity_mag_wall_systole,1));
    stats.vorticity_wall_peak_diastole = max(mean(vorticity_mag_wall_diastole,1));
%  
end