function [out_data, out_name] = Flow2Case(d_path,ser_flow,ser_cine,contour_files,case_file,opts)
% Aaron Hess
% University of Oxford
% last update Feb 2018
% Load flow from dicom series
% perform background correction
% Load cines from a 2D stack
% Load CMR42 contours and match the cine stack
% create flow mask based on the contours
% export to Ensight case format
% opts - structure with options 
% opts.force_reload set true to reload dicom images 
% opts.repeatcycle  set true to repeat cycle for retrogated scans will
%                   double number of frames in case file
% opts.dicom_dir    name of sub folder for dicoms default is 'dicoms'
% NOTE TODO:
% 1. Ensight CASE is in ASCKII, this is slow, must update to binary
% 2. Cine time frames are assumed to be spaced the same as the flow

% [lib_dir,~,~]= fileparts(mfilename('fullpath'));
% lib_dir = [lib_dir,filesep,'/../'];
% addpath([lib_dir,'dcm_orientation_tx']);
% addpath([lib_dir,'cmr42read']);
% addpath([lib_dir,'DicomToolsAH']);


if(~exist('opts','var'))
    opts = [];
end
if(~isfield(opts,'force_reload'))
    opts.force_reload = false;
end
if(~isfield(opts,'repeatcycle'))
    opts.repeatcycle = false;
end
if(~isfield(opts,'dicom_dir'))
    opts.dicom_dir = 'dicoms';
end

rep_cycle = 1;
if(opts.repeatcycle)
    rep_cycle = 2;
end

%% 1. Read dicom files, or load existing mat for flow images


out_data = {};  % nothing to output yet
out_name = {};

mat_file_raw = fullfile(d_path,['flow_raw',num2str(ser_flow),'.mat']);

mat_file_pro = fullfile(d_path,['flow_pro',num2str(ser_flow),'.mat']);
mat_file_cont = fullfile(d_path,['contours',num2str(ser_cine),'.mat']);
mat_file_cine = fullfile(d_path,['flow_cine',num2str(ser_cine),'.mat']);

fprintf('\nLoading flow\n');
if(~exist(mat_file_raw,'file')|| opts.force_reload)
    % Read in data
    fprintf('Importing flow dicoms\n');
    opt = [];
    opt.round_trigger_time = -1;  %auto set to sequence repetition time to gather time frames 
    [im3D, info] = Load3DDicom(fullfile(d_path,opts.dicom_dir),ser_flow,opt);

    % save to matfile
    save(mat_file_raw,'im3D','info','-v7.3');
else
    load(mat_file_raw);
end


% book keeping
im3D = squeeze(im3D);
sz = size(im3D);
NTime = sz(4);
NEncode = sz(5);
flow_dim = 5;
time_dim = 4;

% convert to velocity 
ind = strfind(info.order.contrastsList{NEncode},'_v');
venc  = sscanf(info.order.contrastsList{NEncode}(ind:end),'_v%d')/100; % m/s


flow = squeeze(angle(im3D))*venc/pi();
mag = squeeze(sum(abs(im3D),5));

% Transform to patient (SCT) coordinates, including flow
[~, sl_index] = UpdateList(info.order.sliceList,info.Private_0051_100d);
CRS2SCTflow = calcCRS2SCT(info);  % Note something not quite right here...
CRS2SCTflow(:,4) = CRS2SCTflow*([-1 -1 -sl_index 1])'; % correct for slice which info belongs to, and so that index[1 1 1] is what was [0 0 0]
CRS2SCTflow = CRS2SCTflow/1000; % convert to m (from mm)
CRS2SCTflow(4,4) = 1;

%% 2. reorder flow images for dims to match (1)SAG (R->L)   (2)COR (A->P)  (3)TRA  (F->H)  which is X Y Z in ensight
% NOTE: TO DO, rotate flow vector if slice rotated from these axes

% 'fl3d1_v110in'    'fl3d1_v110ap'    'fl3d1_v110rl'
% also should be rotated... not ure how yet, perhaps look up labels in ice
flow_order = 1:3;
flow_sign = ones(1,1,1,1,3);

NormSCT = cross(info.ImageOrientationPatient(1:3),info.ImageOrientationPatient(4:end));
if((sum(abs(NormSCT)>0)>1)||(sum(abs(info.ImageOrientationPatient(1:3))>0)>1))
    warning('Flow2Case: Image appears rotated, this may not be supported');
end
inPlaneRotation = 0; % not needed for this
[~, ~, orienation] = fGSLCalcPRS(NormSCT, inPlaneRotation);


% for iC = 1:3
%     switch(info.order.contrastsList{iC}(end-1:end))
%         case 'in' %note this is the slice direction
%             flow_order(orienation) = iC;  % orietnation is SAG(1), COR(2), TRA(3)
%             flow_sign(orienation) = 1;
%         case 'fh'
%             flow_order(3) = iC;
%             flow_sign(3) = 1;
%         case 'ap'
%             flow_order(2) = iC;
%             flow_sign(2) = 1;
%         case 'rl'
%             flow_order(1) = iC;
%             flow_sign(1) = 1;
%     end
% end
for iC = 1:3
    switch(info.order.contrastsList{iC}(end-1:end))
        case 'in' %note this is the slice direction
            flow_order(3) = iC;  % orietnation is SAG(1), COR(2), TRA(3)
            if(orienation == 1)  % in prysma data it apears that sagital has a sign reversal for sagital data, must check more
                flow_sign(3) = -1;
            end
        case 'fh'
            flow_order(1) = iC;  % allwas goes first           
            flow_sign(1) = -1;
        case 'ap'
            if(orienation == 3) % TRA
                flow_order(1) = iC;
            else
                flow_order(2) = iC;
            end
            
        case 'rl'
             flow_order(2) = iC;
            
    end
end

flow_in = bsxfun(@times,flow(:,:,:,:,flow_order),flow_sign);

rot_flow = bsxfun(@rdivide,CRS2SCTflow(1:3,1:3),sqrt(sum(CRS2SCTflow(1:3,1:3).^2)));

flow = zeros(size(flow_in));

for i = 1:3
    v = permute(rot_flow(i,:),[1 3 4 5 2]);  % vector in fith dim
    flow(:,:,:,:,i) = sum(bsxfun(@times,flow_in,v),5);
end

%% 3. background phase correct
background_fit = BackgroundCorrection(flow,venc,2);
flow = bsxfun(@minus, flow, background_fit);


%% 4. create a mask using SNR of magnitude
SNR = mean(abs(im3D(:,:,:,:,1)),time_dim)./std(diff(abs(im3D(:,:,:,:,1)),1,time_dim),[],time_dim);
mask = SNR>5;  % SNR Threshold based on magnitude data

flow = bsxfun(@times,flow,mask);


%% 5. anti-aliasing?
% still to do





%% 7. Processed variables
flow_mag = sqrt(sum(flow.^2,flow_dim));

peak_flow = max(flow_mag,[],time_dim);

threshold = 0.1; % m/s
stasis = mask.*sum(flow_mag<threshold,time_dim)/NTime;

% 8. Save results
save(mat_file_pro,'flow','flow_mag','stasis','peak_flow','CRS2SCTflow');

%% 9. Read in cine
iscine = ~isempty(ser_cine);
if(~iscine)
    fprintf('skiping loading cine\n');
else
    fprintf('loading cine\n');
    opt.round_trigger_time = -1;
    if(~exist(mat_file_cine,'file') || opts.force_reload)
        [imCine3D, infoCine] = Load3DDicom(fullfile(d_path,opts.dicom_dir),ser_cine,opt);
        [~, sl_index] = UpdateList(infoCine.order.sliceList,infoCine.Private_0051_100d);
        CRS2SCTCine = calcCRS2SCT(infoCine);
        CRS2SCTCine(:,4) = CRS2SCTCine*([0 0 -sl_index+1 1])'; % correct for slice which info belongs to
        CRS2SCTCine = CRS2SCTCine/1000;  % conver to m (from mm)
        CRS2SCTCine(4,4) = 1;
        save(mat_file_cine, 'imCine3D', 'infoCine', 'CRS2SCTCine');
    else
        load(mat_file_cine)
    end

    [~, sl_index] = UpdateList(infoCine.order.sliceList,infoCine.Private_0051_100d);
    CRS2SCTCine = calcCRS2SCT(infoCine);
    CRS2SCTCine(:,4) = CRS2SCTCine*([0 0 -sl_index+1 1])'; % correct for slice which info belongs to
    CRS2SCTCine = CRS2SCTCine/1000;  % conver to m (from mm)
    CRS2SCTCine(4,4) = 1;

    no_uid = cellfun(@isempty,infoCine.uids);
    infoCine.uids(no_uid) = {'empty'};
    sz_cine = size(imCine3D);
    NPhsCine = size(imCine3D,6);
end
%% 10. Read in controus and generate masks
% match each image to a contour
% NOTE NOT ASSUMING ALL CONTOURS ARE NEEDED
% Must strip all contours not matching our label...
iscontours = (~isempty(contour_files))&&iscine;
if(~( iscontours))
    fprintf('No cine or contours provided\n');
else
    fprintf('loading contours\n');
    if(iscell(contour_files))
        n_cont = length(contour_files);
    else
        contour_files = {contour_files};
        n_cont = 1;
    end

    % if contour file exists load it so we can update it
    if(exist(mat_file_cont,'var') &&~opts.force_reload)
        load(mat_file_cont);
    end

    labels = {};

    for iContFile = 1:n_cont
        [contours, study_uid] = ExtractCMR42Contours(fullfile(d_path,contour_files{iContFile}));

        % add label to this to allow multiple labels
        for iImg = 1:length(contours)
            for iCont = 1:length(contours(iImg).str)
                [labels,ind] = UpdateList(labels,contours(iImg).str(iCont).Label);
            end
        end

        [~,c_nameBase,~] = fileparts(contour_files{iContFile});

        for iLabel = 1:length(labels)
            c_name = [c_nameBase,labels{iLabel}];
        %     cont_match = zeros(length(contours),8);
            cont_vol = zeros(size(infoCine.uids));
            areaToVol = 1e-3*prod(diag(CRS2SCTCine))/4/4;  % units ml (4 as ,Points are x 4 ...)

        %     phsCine2Flow = NTime/NPhsCine;

            BWCine = zeros(size(imCine3D));

            clear('pt_all');
            pt_all{sz_cine(6)} = [];
            dcm_found = false;


            for iImg = 1:length(contours)

                is_match = strfind(infoCine.uids,contours(iImg).uid);
                ind = find(~cellfun(@isempty,is_match));
                % create a label list too...
                if(isempty(ind))
                    warning(['Could not match contours to cine, contour file: ',contour_files{iContFile}]);
                    continue;
                end
                [a, b, c, d, e, f, g] = ind2sub(size(infoCine.uids),ind);  % a is slice d is time frame
        %         cont_match(iImg,2:8) = [a, b, c, d, e, f,g];  
        %         cont_match(iImg,1) = ind;
                for iCont = 1:length(contours(iImg).str)
                    if(isempty(contours(iImg).str(iCont).Label))
                        contours(iImg).str(iCont).Label = 'na';
                    end                    
                    if(strcmp(labels{iLabel},contours(iImg).str(iCont).Label))
                        cont_vol(a,b,c,d,e,f,g) = areaToVol*polyarea(contours(iImg).str(iCont).Points(1,:),contours(iImg).str(iCont).Points(2,:));


                        % convert each controur into 4D flow coordinates then create a mask
                        pt = [contours(iImg).str(iCont).Points(2:-1:1,:)/4; repmat(a,[1,size(contours(iImg).str(iCont).Points,2)]); repmat(1,[1,size(contours(iImg).str(iCont).Points,2)])];
                        %ptFlow = cine2flow*pt;
                        BWCine(:,:,a,b,c,d,e,f,g) = roipoly(imCine3D(:,:,1),pt(2,:),pt(1,:));
                        % now match ECG time stamp?
                        pt_all{d} = [pt_all{d},pt];
                        dcm_found = true;
                    end
                end
            end

            if(~dcm_found)
                warning(['Image match not found for contours labeled: ',labels{iLabel}]);
                continue;
            end
            % remove contrast dimenssion
            if(size(BWCine,8)>1)
                BWCine = sum(BWCine,8);
            end

            vol_t = squeeze(sum(cont_vol(:,1,1,:,1,f),1));

            BWCine = squeeze(BWCine);

            % create a flow mask from contours
            BWFlow = zeros(size(flow_mag));
            [X, Y, Z, T,C] = ind2sub(size(BWFlow),1:numel(BWFlow));

            % Tcine = T/phsCine2Flow;  %this streaches cine contours to fit flow
            Tcine = T;
            Tcine(T>size(BWCine,4))=size(BWCine,4); % this assumes same timing just shorter cycle
            inCine = CRS2SCTCine\CRS2SCTflow*[X;Y;Z;C];
            inCine(1,inCine(1,:)<1) = 1;
            inCine(2,inCine(2,:)<1) = 1;
            inCine(3,inCine(3,:)<1) = 1;
            inCine(1,inCine(1,:)>size(BWCine,1)) = size(BWCine,1);
            inCine(2,inCine(2,:)>size(BWCine,2)) = size(BWCine,2);
            inCine(3,inCine(3,:)>size(BWCine,3)) = size(BWCine,3);

            inCine = round(inCine);
            BWFlow(1:end) = BWCine(sub2ind(size(BWCine),inCine(1,:),inCine(2,:),inCine(3,:),round(Tcine)));

            cont.(c_name).BWCine = BWCine;
            cont.(c_name).BWFlow = BWFlow;
            cont.(c_name).pt_all = pt_all;

            cont.(c_name).vol_t = vol_t;


            var.([c_name,'Mask']) = repmat(BWFlow,[1 1 1 rep_cycle]);
        end

    end
    save(mat_file_cont,'cont','CRS2SCTCine');
end  % end is contours

%% Write to Ensight case format

timesteps = info.order.phsList/1000; % 
if(opts.repeatcycle)
    timesteps = [timesteps,timesteps+timesteps(end)];
end
var.stasis = stasis;
var.peakFlow = peak_flow;  % make 1D
var.velocity = repmat(flow,[1 1 1 rep_cycle 1]);
var.PCMRIMag = repmat(mag ,[1 1 1 rep_cycle 1]);



% Match cine time frames to flow
% copy last frame to end
if(iscine)
    sz_cine = size(squeeze(imCine3D));
    nt = sz_cine(4);
    var_cine.cine = zeros([sz_cine(1:3),sz(4)*rep_cycle]);
    var_cine.cine(:,:,:,1:nt) = squeeze(imCine3D(:,:,:,1,1,:,1,2));
    for iT=nt+1:sz(4)*rep_cycle
        ind = mod(iT-1,sz(4))+1;
        if(ind>nt) 
            ind = nt; 
        end
        var_cine.cine(:,:,:,iT) = squeeze(imCine3D(:,:,:,1,1,ind,1,2));
    end
else
    var_cine = [];
    CRS2SCTCine = [];
end
if(~iscontours)
    cont = [];
end

cfile = [d_path,filesep,case_file,filesep,case_file,'.case'];
[dirCase,~,~] = fileparts(cfile);
mkdir(dirCase);
fprintf('Exporting to case file: %s\n',cfile);
Ensight4DFlow(cfile,timesteps,CRS2SCTflow,var,CRS2SCTCine,var_cine,CRS2SCTCine,cont);


