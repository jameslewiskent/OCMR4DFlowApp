% ATH Sep 2019, check that all flow datasets have 30 time frames
% Check the following in the data
% 1. all flow has 30 time frames, this garantees they all were temporally
% binned correctly


DATA_LIST_GPA;

n_studies = length(dat);
evaluate_studies = {''}; % run this to process this list of studies, seperate with commas
 evaluate_studies = []; % run this to process all studies

if (isempty(evaluate_studies))
    run_studies  =1:n_studies;
else
    run_studies = 1:length(evaluate_studies);
end

for i_proc = 1:run_studies(end)
    if(isempty(evaluate_studies))
        i_study = i_proc;
    else
        for i_test = 1:n_studies
            
            if(~isempty(dat(i_test))&&~isempty((strfind(dat(i_test).d_path,evaluate_studies{i_proc}))))
                i_study = i_test;
                break;
            elseif(i_test==n_studies)
                i_study = [];  % could not find match
                fprintf('\nNo Match found for study %s\n',evaluate_studies{i_proc});
            end
        end
    end
    if(isempty(i_study)||isempty(dat(i_study)))
        continue; %ignore this study as no match
    end

 
    fprintf('\nProcessing study %s\n',dat(i_study).d_path);
     
    d_path = dat(i_study).d_path;
    contour_file = dat(i_study).contour_file;
    ser_2Dcines = dat(i_study).ser_2Dcines;
    ser_flow = dat(i_study).ser_flow;
    
    mat_file_pro = fullfile(d_path,['SCAN\flow_pro',num2str(ser_flow),'.mat']);
    if(exist(mat_file_pro,'file'))
        mf = matfile(mat_file_pro);
        nt_flow(i_proc) =size(mf,'flow',4);

        if(nt_flow(i_proc)>30)
            fprintf('Study %s not consistant\n',dat(i_study).d_path);
        end
    else
        fprintf('Study not found %s\n',dat(i_study).d_path);
    end
end