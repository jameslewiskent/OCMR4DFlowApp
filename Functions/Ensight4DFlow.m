function Ensight4DFlow(fname,timesteps,CRS2SCTflow,vars,CRS2SCTCine,var_cine,CRS2SCTContours,contours)
    n_contours_to_write = 10;
    %[lib_dir,~,~]= fileparts(mfilename('fullpath'));
    %addpath([lib_dir,filesep,'EnsightWriter']);
    
    %ensure even spacing of timesteps
    dT =  round(1000*(timesteps(end)-timesteps(1))/(length(timesteps)-1))/1000;
    timesteps = (0:length(timesteps)-1)*dT;
    cfile = CaseFile(timesteps);
    var_names = fieldnames(vars);
    sz = size(vars.(var_names{1}));
    flowpart = cfile.CreatePart('flow4D');
    flowpart.CreateBlockGeometry(CRS2SCTflow,sz(1:3));
    
    
    for iV = 1:length(var_names)
        flowpart.SetVariableTime(var_names{iV},vars.(var_names{iV}));
    end
    
    % create a cine part... must have same number time frames
    if(~isempty(CRS2SCTCine))
        cinepart = cfile.CreatePart('cine3D');
        sz_cine = size(var_cine.cine);
        cinepart.CreateBlockGeometry(CRS2SCTCine,sz_cine(1:3));
        var_names = fieldnames(var_cine);
        for iV = 1:length(var_names)
            cinepart.SetVariableTime(var_names{iV},var_cine.(var_names{iV}));
        end
    end
    %% write contours into their own dataset based on cine CRS2SCR
    if(~isempty(contours))
        cont_names = fieldnames(contours);
        if(n_contours_to_write>length(cont_names))
            n_contours_to_write = length(cont_names);
        end
        for iCont = 1:n_contours_to_write
            contpart = cfile.CreatePart(cont_names{iCont});

            % Note the times of the contours are differnt from flow, must resolve
            frm = 1;
            for iT=1:sz(4)
                if(length(contours.(cont_names{iCont}).pt_all)<iT) || isempty(contours.(cont_names{iCont}).pt_all{iT})
        %             frm = last_frm;
                else
                    frm = iT;
                end
                vert_cont = (CRS2SCTContours*contours.(cont_names{iCont}).pt_all{frm})';
                contpart.SetGeometry(iT,vert_cont(:,1),vert_cont(:,2),vert_cont(:,3));
                TRI = delaunay(vert_cont(:,1),vert_cont(:,2),vert_cont(:,3));
                contpart.SetConnectivity(iT,TRI,'tetra4');
            end
        end
    end
    cfile.WriteFile(fname);
%     
end