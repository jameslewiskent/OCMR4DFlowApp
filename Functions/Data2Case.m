function Data2Case(case_dir,case_name,parts,app)
% Aaron Hess
% University of Oxford
% parts is a structer list of all parts 
% parts(i).vars, parts(i).timesteps, parts(i).CRS2SCT, parts(i).name
% time steps are synchronised to the first part
% April 2018

    %[lib_dir,~,~]= fileparts(mfilename('fullpath'));
    %addpath([lib_dir,filesep,'EnsightWriter']);
    
    %ensure even spacing of timesteps
    dT =  round(1000*(parts(1).timesteps(end)-parts(1).timesteps(1))/(length(parts(1).timesteps)-1))/1000;
    timesteps = (0:length(parts(1).timesteps)-1)*dT;
    
    cfile = CaseFile(timesteps);
    
    
    for i_part = 1:length(parts)
        vars = parts(i_part).vars;
        var_names = fieldnames(vars);
        sz = size(vars.(var_names{1}));
        sz(length(sz)+1:4) = 1;

        part = cfile.CreatePart( parts(i_part).name);
        part.CreateBlockGeometry(parts(i_part).CRS2SCT,sz(1:3));

        for iV = 1:length(var_names)
            nT = size(vars.(var_names{iV}),4);
            if (nT>1)&&(nT~=length(timesteps))
                % map timesteps to timesteps of first part
                dTpart = max(parts(i_part).timesteps)/length(parts(i_part).timesteps);
                iv_ind  = round(timesteps/dTpart);
                iv_ind(iv_ind<1) = 1;
%                 iv_ind(iv_ind>length(timesteps)) = length(timesteps);
                iv_ind(iv_ind>nT) = nT;
                var = vars.(var_names{iV})(:,:,:,iv_ind);
                part.SetVariableTime(var_names{iV},var);
            else
                part.SetVariableTime(var_names{iV},vars.(var_names{iV}));
            end
        end
    end
    
    cfile_name = [case_dir,filesep,case_name,'.case'];
    mkdir(case_dir);
    
    cfile.WriteFile(cfile_name,app);
end