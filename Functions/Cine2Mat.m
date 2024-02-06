function [out_data, out_file] = Cine2Mat(d_path,ser_cine,opts)
% Aaron Hess
% University of Oxford
% April 2018
% Read in a 2D/3D dicom and create 
% Load cines from a 2D stack in ser_cine 

    if(~exist('opts','var'))
        opts = [];
    end
    if(~isfield(opts,'force_reload'))
        opts.force_reload = false;
    end
    if(~isfield(opts,'round_trigger_time'))
        opts.round_trigger_time = -1;
    end
    mat_file_cine = fullfile(d_path,['flow_cine',num2str(ser_cine),'.mat']);


    fprintf('loading cine\n');
    opt.round_trigger_time = opts.round_trigger_time;
    if(~exist(mat_file_cine,'file') || opts.force_reload)
        [imCine, infoCine] = Load3DDicom(fullfile(d_path,opts.dicom_dir),ser_cine,opt);
        [~, sl_index] = UpdateList(infoCine.order.sliceList,infoCine.Private_0051_100d);
        CRS2SCTCine = calcCRS2SCT(infoCine);
        CRS2SCTCine(:,4) = CRS2SCTCine*([0 0 -sl_index+1 1])'; % correct for slice which info belongs to
        CRS2SCTCine = CRS2SCTCine/1000;  % conver to m (from mm)
        CRS2SCTCine(4,4) = 1;
        save(mat_file_cine, 'imCine', 'infoCine', 'CRS2SCTCine');
    else
        load(mat_file_cine)
    end

%     [~, sl_index] = UpdateList(infoCine.order.sliceList,infoCine.Private_0051_100d);
%     CRS2SCTCine = calcCRS2SCT(infoCine);
%     CRS2SCTCine(:,4) = CRS2SCTCine*([0 0 -sl_index+1 1])'; % correct for slice which info belongs to
%     CRS2SCTCine = CRS2SCTCine/1000;  % conver to m (from mm)
%     CRS2SCTCine(4,4) = 1;

%     no_uid = cellfun(@isempty,infoCine.uids);
%     infoCine.uids(no_uid) = {'empty'};
%     sz_cine = size(imCine);
%     NPhsCine = size(imCine,6);
    
    
    out_data.imCine = imCine;
    out_data.infoCine = infoCine;
    out_data.CRS2SCTCine = CRS2SCTCine;
    
    out_file = mat_file_cine;
end