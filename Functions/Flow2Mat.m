function [out_data, out_name, info] = Flow2Mat(d_path,ser_flow,opts,app,progbar)
% Aaron Hess
% University of Oxford
% last update April 2018
% Load flow from dicom series

% perform background correction
% opts - structure with options 
% opts.force_reload set true to reload dicom images 
% opts.repeatcycle  set true to repeat cycle for retrogated scans will
%                   double number of frames in case file
% opts.dicom_dir    name of sub folder for dicoms default is 'dicoms'
% opts.phaseunwrapping  1 = on, 0 = off  (default = 1)
% NOTE TODO:
% 1. Ensight CASE is in ASCKII, this is slow, must update to binary
% 2. Cine time frames are assumed to be spaced the same as the flow

%%
% Load CMR42 contours and match the cine stack
% create flow mask based on the contours
% export to Ensight case format

%[lib_dir,~,~]= fileparts(mfilename('fullpath'));
%lib_dir = [lib_dir,filesep,'/../'];
%addpath([lib_dir,'dcm_orientation_tx']);
%addpath([lib_dir,'cmr42read']);
%addpath([lib_dir,'DicomToolsAH']);
%addpath([lib_dir,'DivFreeWavelet/Func_DFW']);

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
    opts.dicom_dir = 'dicom';
end

if(~isfield(opts,'phaseunwrapping'))
    opts.phaseunwrapping = 1; 
end

if(~isfield(opts,'noisereduction'))
    opts.noisereduction = 0;
end

% rep_cycle = 1;
% if(opts.repeatcycle)
%     rep_cycle = 2;
% end

% JK added for dealing with contours on phase image
if(~isfield(opts,'CountoursonPhaseImage'))
    opts.CountoursonPhaseImage = 0;
end
%% 1. Read dicom files, or load existing mat for flow images


out_data = {};  % nothing to output yet

mat_file_raw = fullfile(d_path,['flow_raw',num2str(ser_flow),'.mat']);

mat_file_pro = fullfile(d_path,['flow_pro',num2str(ser_flow),'.mat']);

%mat_file_cine = fullfile(d_path,['flow_cine',num2str(ser_cine),'.mat']);

fprintf('\nLoading flow\n');
if(~exist(mat_file_raw,'file')|| opts.force_reload)
    % Read in data
    fprintf('Importing flow dicom\n');
    opt = opts;
    opt.round_trigger_time = -1;  %auto set to sequence repetition time to gather time frames 
    
    [im3D, info] = Load3DDicom(fullfile(d_path,opts.dicom_dir),ser_flow,opt,app);

            if opts.CountoursonPhaseImage~=0
                % JK added. Contours performed on phase image, replace
                % magnitude info.uids with those from contoured phase image
                % This is a bit of a WIP hack
                ser_flow(1) = ser_flow(1) - 1;
                [~, tempinfo] = Load3DDicom(fullfile(d_path,opts.dicom_dir),ser_flow,opt,app);
                info.uids = tempinfo.uids;
            end
    
    % save to matfile
    save(mat_file_raw,'im3D','info','-v7.3');
else
    load(mat_file_raw,'im3D','info'); % JK added ,'im3D','info'
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
SliceL = uint8(info.Private_0051_100d);
[~, sl_index] = UpdateList(info.order.sliceList,SliceL); %info.Private_0051_100d
CRS2SCTflow = calcCRS2SCT(info);  % Note something not quite right here...
CRS2SCTflow(:,4) = CRS2SCTflow*([-1 -1 -sl_index 1])'; % correct for slice which info belongs to, and so that index[1 1 1] is what was [0 0 0]
CRS2SCTflow = CRS2SCTflow/1000; % convert to m (from mm)
CRS2SCTflow(4,4) = 1;

%% 2. reorder flow images for dims to match (1)SAG (R->L)   (2)COR (A->P)  (3)TRA  (F->H)  which is X Y Z in ensight
% NOTE: TO DO, rotate flow vector if slice rotated from these axes

% 'fl3d1_v110in'    'fl3d1_v110ap'    'fl3d1_v110rl'
% also should be rotated... not sure how yet, perhaps look up labels in ice
flow_order = 1:3;
flow_sign = ones(1,1,1,1,3);

NormSCT = cross(info.ImageOrientationPatient(1:3),info.ImageOrientationPatient(4:end));
if((sum(abs(NormSCT)>0)>1)||(sum(abs(info.ImageOrientationPatient(1:3))>0)>1))
     uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'Image appears rotated, this may not be supported.'.','An error occured');
     app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Image appears rotated, this may not be supported. '); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
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
for iC = 1:length(info.order.contrastsList)
    switch(info.order.contrastsList{iC}(end-1:end))
        case 'in' % note this is the slice direction
            flow_order(3) = iC;  % orientation is SAG(1), COR(2), TRA(3)
            if(orienation == 1)  % in prisma data it apears that sagital has a sign reversal for sagital data, must check more
                flow_sign(3) = -1;
            end
        case 'fh'
            flow_order(1) = iC;  % always goes first           
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
    
                if progbar.CancelRequested
                    break
                end
end

flow_in = bsxfun(@times,flow(:,:,:,:,flow_order),flow_sign);

rot_flow = bsxfun(@rdivide,CRS2SCTflow(1:3,1:3),sqrt(sum(CRS2SCTflow(1:3,1:3).^2)));

flow = zeros(size(flow_in));

for i = 1:3
    v = permute(rot_flow(i,:),[1 3 4 5 2]);  % vector in fifth dim
    flow(:,:,:,:,i) = sum(bsxfun(@times,flow_in,v),5);
end

%% 3. anti-aliasing?
% still to do
if(opts.phaseunwrapping)
    %addpath('unwrap')
    for i =  1:3  % for each direction, unwrap
        n_u4 = unwrap_4D(flow(:,:,:,:,i)/(venc/pi()));
        flow(:,:,:,:,i) = flow(:,:,:,:,i) + 2*pi .* double(n_u4)*(venc/pi());
    end
end

%% 4. background phase correct
background_fit = BackgroundCorrection(flow,venc,2);
flow = bsxfun(@minus, flow, background_fit);


%% 5. create a mask using SNR of magnitude
SNR = mean(abs(im3D(:,:,:,:,1)),time_dim)./std(diff(abs(im3D(:,:,:,:,1)),1,time_dim),[],time_dim);

dT = info.order.phsList(2) - info.order.phsList(1);
if(dT<40)
    mask = SNR>5;  % SNR Threshold based on magnitude data
else
    mask = SNR>(5*25/dT);
end
flow_unmasked = flow;
flow = bsxfun(@times,flow,mask);





%% 6. denoising
if(opts.noisereduction)
    spins = 2;              % Number of cycle spinning per dimension
    isRandShift = 1;        % Use random shift
    minSize = 8*ones(1,3);  % Smallest wavelet level size
    res = [1,1,1];
%     flow_dfw = zeros(size(flow));

    for it = 1:size(flow,4)
        [vxDFW_sms,vyDFW_sms,vzDFW_sms] = dfwavelet_thresh_SURE_MAD_spin(flow(:,:,:,it,1),flow(:,:,:,it,2),flow(:,:,:,it,3),minSize,res,spins,isRandShift);
        flow_unmasked(:,:,:,it,1) = vxDFW_sms;
        flow_unmasked(:,:,:,it,2) = vyDFW_sms;
        flow_unmasked(:,:,:,it,3) = vzDFW_sms;
    end

end

% Vorticity
ps  = info.PixelSpacing/1000;

[X, Y, Z] = meshgrid((1:size(flow,2))*ps(1),(1:size(flow,1))*ps(2),(1:size(flow,3))*info.SliceThickness/1000);
curlx = zeros(size(flow,1),size(flow,2),size(flow,3),size(flow,4));
curly = zeros(size(flow,1),size(flow,2),size(flow,3),size(flow,4));
curlz = zeros(size(flow,1),size(flow,2),size(flow,3),size(flow,4));
cav = zeros(size(flow,1),size(flow,2),size(flow,3),size(flow,4));
for it = 1:size(flow,4)
    [curlx(:,:,:,it),curly(:,:,:,it),curlz(:,:,:,it),cav(:,:,:,it)] = curl(X,Y,Z,flow_unmasked(:,:,:,it,1),flow_unmasked(:,:,:,it,2),flow_unmasked(:,:,:,it,3)) ;
end

vorticity = sqrt(curlx.^2+curly.^2+curlz.^2);
vorticity_map = mask.*sum(vorticity,time_dim)/NTime;

%% lambda2 for vorticity

lambda2 = Lambda2Calc(flow_unmasked,[ps([2,1]);info.SliceThickness/1000]);

%% 7. Processed variables
flow_mag = sqrt(sum(flow.^2,flow_dim));
pcmri= mag.*(sum(flow_unmasked.^2,flow_dim).^(0.2));
peak_flow = max(flow_mag,[],time_dim);
threshold = 0.1; % m/s size*
stasis = mask.*sum(flow_mag<threshold,time_dim)/NTime;
% % flow_mag_systole=flow_mag(:,:,:,[1:26]);
% % flow_mag_diastole=flow_mag(:,:,:,[27:30]);
% % vorticity_systole = vorticity(:,:,:,[1:26]);
% % vorticity_diastole = vorticity(:,:,:,[27:30]);
% 
% %time analysis
% stasis_systole = mask.*sum(flow_mag_systole<threshold,time_dim)/26;
% stasis_diastole = mask.*sum(flow_mag_diastole<threshold,time_dim)/4
% peak_flow_systole = max(flow_mag_systole,[],time_dim);
% peak_flow_diastole = max(flow_mag_diastole,[],time_dim);

   
% 8. Save results
save(mat_file_pro,'flow','flow_mag','stasis','peak_flow','CRS2SCTflow','info','mag','pcmri','flow_unmasked','vorticity','vorticity_map','cav','lambda2');
out_data.flow = flow;
out_data.flow_mag = flow_mag;
% out_data.flow_mag_systole = flow_mag_systole;
% out_data.flow_mag_diastole = flow_mag_diastole;
out_data.stasis = stasis;
% out_data.stasis_systole = stasis_systole;
% out_data.stasis_diastole = stasis_diastole;
% out_data.peak_flow_systole = peak_flow_systole;
% out_data.peak_flow_diastole = peak_flow_diastole;
out_data.peak_flow = peak_flow;
out_data.CRS2SCTflow = CRS2SCTflow;
out_data.mag = mag;
out_data.info = info;
out_data.pcmri = pcmri;
out_data.flow_unmasked = flow_unmasked;
out_data.vorticity = vorticity;
% out_data.vorticity_systole = vorticity_systole;
% out_data.vorticity_diastole = vorticity_diastole;
out_data.vorticity_map = vorticity_map;
out_data.cav = cav;
out_data.lambda2 = lambda2;

out_name = mat_file_pro;


