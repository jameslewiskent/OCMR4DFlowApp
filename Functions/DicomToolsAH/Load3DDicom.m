function [im3D, info] = Load3DDicom(data_dir,data_series,opt,app)
%  Load3DDicom(data_dir,data_series,opt)
% Sort dicom directory according to measured channel, echo, and slice no
% Magnitude and phase are combined into a complex output
% output structure is NX NY NSlc NCh NTE
% options struct opt.
%   ignore_trigger_time (bool) - disable trigger time detection
%   round_trigger_time (int) - round trigger time in ms
% Dimensions are
%     Read Phase Slice Ch TE Repetition Contrast
% valid options (opt) are
% ignore_trigger_time  - do not care about ECG
% round_trigger_time  - round to nearst ms .eg for 30ms resolution = 30 ms
%                     - if -1 then set it to TR
% ignore_duplicates - if two files map to the same 3D data, ignore the
%                     second of them

progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Searching DICOM files'); pause(0.1)

NDim = 6;

if(exist('opt','var'))
    if(~isfield(opt,'ignore_trigger_time'))
        opt.ignore_trigger_time = 0;
    end
%     if(~isfield(opt,'round_trigger_time'))
%         opt.round_trigger_time = 1;
%     end
    if(~isfield(opt,'ignore_duplicates'))
        opt.ignore_duplicates = 0;
    end
    
else
    opt.ignore_trigger_time = 0;
%     opt.round_trigger_time = 1;
    opt.ignore_duplicates = 0;
end
fprintf('Load3DDicom settings:\nopt.ignore_trigger_time = %g\nopt.round_trigger_time=%g\nopt.ignore_duplicates=%g\n',opt.ignore_trigger_time,opt.round_trigger_time,opt.ignore_duplicates);
%
%fprintf('\n1/3 Searching for Dicom files');
[dicomData] = processDicomDirRecursive(app,data_dir, '*'); 
    
NSer = length(dicomData.study(1).series);
name_list = cell(NSer,1);

for iSer = 1:NSer
    s = dicomData.study(1).series(iSer);
    name_list{iSer} = [num2str(s.SeriesNumber), '_', s.SeriesDescription,'_',num2str(length(s.instance))];    
end
if (~exist('data_series','var') || isempty(data_series))
    % Match series descriptions instread of selection list (future feature) in setstruct(1,4).SeriesDescription
    [ser_ind,~] = listdlg('PromptString','Select mag and phase uncombined 3D datasets:',...
                        'SelectionMode','Multiple',...
                        'ListString',name_list);
else
    ser_ind = zeros(size(data_series));
     for i = 1:length(data_series)
         ser_ind(i) = find([dicomData.study(1).series.SeriesNumber] == data_series(i));
     end
end
% Check the two series are the same length
N_IMG = length(dicomData.study(1).series(ser_ind(1)).instance);
for iSer = 2:length(ser_ind)
    if ( N_IMG ~= length(dicomData.study(1).series(ser_ind(iSer)).instance))
        fprintf('\nSeries lengths are inconsistant, for ser=%g found %g, expecting %g\n',ser_ind(iSer),length(dicomData.study(1).series(ser_ind(iSer)).instance),N_IMG);
        error('Failed DICOM data load');
    end
end

% Determine size of dataset
% info.size = [NRead NPhase NSlice NCh NTE NContrast];
info.size = zeros(5,1);

N_SER = length(ser_ind);
TOTAL_IMGS = N_IMG*N_SER;

ch_id = zeros(TOTAL_IMGS,1);
echo_no = zeros(TOTAL_IMGS,1);
slice_no = zeros(TOTAL_IMGS,1);
rep_no =  zeros(TOTAL_IMGS,1);
uids = cell(TOTAL_IMGS,1);
cont_no = zeros(TOTAL_IMGS,1);
phs_no = zeros(TOTAL_IMGS,1);
is_phs = zeros(TOTAL_IMGS,1);

triger_time = zeros(TOTAL_IMGS,1);

im2d_in = dicomread(dicomData.study(1).series(ser_ind(1)).instance(1).Filename);
info = dicominfo(dicomData.study(1).series(ser_ind(1)).instance(1).Filename);
opt.shift_cycle = 0;

% ATH 02/2018 auto setting of grouping to sequence repetition time
% if(opt.round_trigger_time==-1)
%     if(info.RepetitionTime<info.NominalInterval/info.CardiacNumberOfImages)
%         opt.round_trigger_time = info.RepetitionTime;
%         opt.shift_cycle = -30; % prospective timeings are offset by navigator delay
%     else
%         if(info.NominalInterval>100)
%             opt.round_trigger_time = info.NominalInterval/info.CardiacNumberOfImages;
%         else
%             opt.round_trigger_time = 1;
%         end
%         
%     end
% end
    
info.size(1:2) = size(im2d_in);

img_in = zeros([info.size(1:2) TOTAL_IMGS]);  % read x phase x imgs x magPhs
teList = [];
chList = {};
sliceList = {};
repList = [];
%phaseList = [];  % cardiac phase either on TriggerTime or?? not sure
%setList = {};    % encoding set, not sure how to determin this yet....
phsList = []; 
contList = {};

%fprintf('\r2/3 Sorting dicom files');
%progressbar(0);
progbar.Value = 1;
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Sorting DICOM files'); pause(0.1)



for  iImg = 1:TOTAL_IMGS  % loop over mag and phase if there are both
    
    iInst = mod(iImg,N_IMG)+1;
    iSer = ceil(iImg/N_IMG);
    s = dicomData.study(1).series(ser_ind(iSer));
        
    img_in(:,:,iImg) = dicomread(s.instance(iInst).Filename);
    info = dicominfo(s.instance(iInst).Filename);

    if(~isfield(info,'TriggerTime'))
        opt.ignore_trigger_time = 1;
    end

    is_phs(iImg) = ((~isempty(strfind(info.ImageType,'\P\'))) || strcmp(info.ImageType(end-1:end),'\P'));
    
    % [out_list, list_ind] = UpdateList(in_list, var)
    try
    [chList, ch_index] = UpdateList(chList,info.Private_0051_100f);
    [sliceList, sl_index] = UpdateList(sliceList,info.Private_0051_100d);
    catch
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - WARNING! Failed in Load3DDicom.m. Trying again.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    % Replaced Load3DDicom_DK.m 
    [chList, ch_index] = UpdateList(chList,uint8(info.Private_0051_100d));
    [sliceList, sl_index] = UpdateList(sliceList,uint8(info.Private_0021_114f));
    end
    
    [teList, te_index] = UpdateList(teList,info.EchoTime);
    [repList, rep_index] = UpdateList(repList,info.AcquisitionNumber); 
    [contList, cont_index] = UpdateList(contList,info.SequenceName); 
%     if(~opt.ignore_trigger_time)
%         [phsList, phs_index] = UpdateList(phsList,opt.round_trigger_time*round((info.TriggerTime+opt.shift_cycle)/opt.round_trigger_time)); 
%     else
%         [phsList, phs_index] = UpdateList(phsList,0); 
%     end
%     if(info.TriggerTime>max_trigger_found)
%         max_trigger_found = info.TriggerTime;
%     end

    ch_id(iImg) = ch_index;   % info.Private_0051_100f
    echo_no(iImg) = te_index;  % info.EchoTime
    slice_no(iImg) = sl_index; % info.Private_0051_100d
    rep_no(iImg) = rep_index;
    phs_no(iImg) = 0;%phs_index; % now done below
    cont_no(iImg) = cont_index;
    uids{iImg}   = info.SOPInstanceUID;
    triger_time(iImg) = info.TriggerTime; %TODO USE THIS TO BIN TIME FRAMES RETROSPECTIVLY!!!!
    
    if(mod(iImg,10)==0)
        progbar.Value = (iImg/TOTAL_IMGS);
    end
    
end

% if(abs(max_trigger_found-max(phsList)-info.NominalInterval)>opt.round_trigger_time)
%     fprintf('\nWarning: cardiac cycle possibly incorectly binned\nNumber time points found = %g\nmaximum time frame = %g\nDicom reported nominal interval = %g\n',length(phsLst),max_trigger_found,info.NominalInterval);
%     fprintf('Rebinning again now\n')
if(~opt.ignore_trigger_time)
    triger_time = triger_time - min(triger_time);  % correct any offsets in prospective
    round_trigger_time = (max(triger_time))/(info.CardiacNumberOfImages-1);
    phsList2 = [];
    for iImg=1:TOTAL_IMGS
        [phsList2, phs_index2] = UpdateList(phsList2,round_trigger_time*round((triger_time(iImg))/round_trigger_time)); 
        phs_no(iImg) = phs_index2;
    end
end
    
NSlice = length(sliceList);
NCh = length(chList);
NTE = length(teList);
NREP = length(repList);
NPHS = length(phsList2);
NCont = length(contList);


sliceListNew = zeros(1,length(sliceList));
% Convert slice list into number
for iSlc = 1:NSlice
    str = sliceList{iSlc}(4:end);
    str = strrep(str,'F','-');
    str = strrep(str,'H','+');    
    str = strrep(str,'L','-');
    str = strrep(str,'R','+');
    % not sure of sign of below
    str = strrep(str,'A','-');
    str = strrep(str,'P','+');
    
    sliceListNew(iSlc) = str2double(str);
end
[~, slro] = sort(sliceListNew);
[~,slroI] = sort(slro); % invert sorting index for reordering

%FIX PHASE LIST ORDER
[~, phro] = sort(phsList2);
[~,phroI] = sort(phro); % invert sorting index for reordering

info.size(1:2) = size(im2d_in);
info.size([3,4,5,6,7,8]) = [NSlice NCh NTE NPHS NREP NCont];
info.TE = teList;
info.order.sliceList = sliceList(slro);
info.order.chList = chList;
info.order.contrastsList = contList;
info.order.phsList=phsList2;
info.order.teList=teList;

phase_dim = 1 + (sum(is_phs)>0);
info.uids = cell([info.size(3:end), phase_dim]);

progbar.Value = 1;
%fprintf('\r3/3 Collating dicom files');
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Collating DICOM files'); pause(0.1)


%
im3D = ones(info.size);
cnt3DM= zeros(info.size(3:end));
cnt3DP= zeros(info.size(3:end));
% [slind chind teind] = meshgrid(1:NSlice, 1:NCh, 1:NTE);
% im3D(:,:,slice_no(:,1),ch_id(:,1),echo_no(:,1)) = img_in(:,:,:,1);
% Just loop it, its faster
for iImg = 1:TOTAL_IMGS
    %im3D(:,:,slice_no(iImg),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg)) = img_in(:,:,iImg).*exp(2i*pi*(img_in(:,:,iImg)-2048)/4096);
    if(opt.ignore_duplicates )
        if(is_phs(iImg))
            if(~cnt3DP(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)))
                im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)).*exp(2i*pi*(img_in(:,:,iImg)-2048)/4096);
                cnt3DP(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = cnt3DP(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg))+1;
                info.uids(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg),is_phs(iImg)+1) = uids(iImg);
            end
        else
            if(~cnt3DM(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)))        
                im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)).*img_in(:,:,iImg);
                cnt3DM(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = cnt3DM(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) +1;
                info.uids(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg),is_phs(iImg)+1) = uids(iImg);
            end
        end
    else
        if(is_phs(iImg))
            im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)).*exp(2i*pi*(img_in(:,:,iImg)-2048)/4096);
            cnt3DP(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = cnt3DP(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg))+1;
        else            
            im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = im3D(:,:,slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)).*img_in(:,:,iImg);
            cnt3DM(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) = cnt3DM(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg)) +1;
        end
        info.uids(slroI(slice_no(iImg)),ch_id(iImg),echo_no(iImg),phs_no(iImg),rep_no(iImg),cont_no(iImg),is_phs(iImg)+1) = uids(iImg);
    end
    
    
    if(mod(iImg,10)==0)
        progbar.Value = (iImg/TOTAL_IMGS);
    end
end

progbar.Value = 1;
if(sum(cnt3DM(:)>1)>0) || (sum(cnt3DP(:)>1)>0)
    a = 1; % allow debug bypass
    if(a)
        uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'More than one image for same index found, likely a problem with the data!','An error occured');
        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - More than one image for same index found, likely a problem with the data.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
        error('Load3DDicom.m: More than one image for same index found, likely a problem with the data!');
    end
end
end



