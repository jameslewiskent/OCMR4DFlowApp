function ConvertCineStack_To3D(dir_in,segment_file,reference_slice)
% Aaron Hess
% Universtiy of Oxford
% Septermber 2014
% Generate a 3D dicom series based on SA cine stack
% and a second series based on a segment mask
% Segment mask is further aligned from slice to slice to ensure a
% contiguous volume, a new segment file is written as
% old_name_corrected.mat
% Usage ConvertCineStack_To3D(dir_cine,segment_file,reference_slice)
% segment_file = 'G:\Data\4DFlow_WholeHeart\DCMIV005\DCMIV005_valign.mat';
% dir_in = 'G:\Data\4DFlow_WholeHeart\DCMIV005\CINES';
% reference_slice (optional) if not specify, it will automatically be chosen

if ~exist('reference_slice','var')
    reference_slice = -1;
end

dir_out = [dir_in, filesep,'../my_files/'];
mkdir(dir_out);

[path name ext] = fileparts(segment_file);
segment_out = [path filesep name '_corrected' ext];

[mfilepath,NAME,EXT] = fileparts(mfilename('fullpath')) ;

addpath([mfilepath, filesep, 'Rodgers_dcm_read']);

[dicomData] = processDicomDirRecursive(dir_in, '*'); 
    
NSer = length(dicomData.study(1).series);
name_list = cell(NSer,1);

for iSer = 1:NSer
    s = dicomData.study(1).series(iSer);
    name_list{iSer} = [num2str(s.SeriesNumber), '_', s.SeriesDescription,'_',num2str(length(s.instance))];
end

%% Match series descriptions instread of selection list (future featuer) in setstruct(1,4).SeriesDescription
[ser_ind,v] = listdlg('PromptString','Select all SA series:',...
                    'SelectionMode','Multiple',...
                    'ListString',name_list);
          
%% Load segment file
load(segment_file);  % Contains im info preview and setstruct


%%
% 
% Convert stack of 2D to 3D dicom series

NZ_SA = length(ser_ind); % Number of images in stack
dcm_infSA = cell(1,NZ_SA);

for i = 1:NZ_SA
    % Need logic to sort these along Z and match to segment info
    dcm_infSA{i} = dicominfo(dicomData.study(1).series(ser_ind(i)).instance(1).Filename);
end

%% Match dicom series to segment series
% Chech that all SA dcm_inf have same orientation
IS_ZERO = 1e-4;
orientation_1 = dcm_infSA{1}.ImageOrientationPatient;
isSame = true;
position_z = zeros(1,NZ_SA);
for i = 1:NZ_SA
    if(sum(abs(orientation_1 - dcm_infSA{i}.ImageOrientationPatient)) > IS_ZERO)
        isSame = false;
    end
    position_z(i) = ...
        cross(dcm_infSA{i}.ImageOrientationPatient(1:3),dcm_infSA{i}.ImageOrientationPatient(4:end))'*dcm_infSA{i}.ImagePositionPatient;
end

% Ensure all are unique
[pos_z_out ind_in ind_out] = unique(position_z,'last');

if length(pos_z_out) ~= length(position_z)
    warning('ConvertMaskToStack: repeated SA dicoms found, latest series used');
end

if (mean(diff(pos_z_out)) - (pos_z_out(2) - pos_z_out(1))) > IS_ZERO
    error('ConvertMaskToStack: Missing slice or uneven slice spacing')
end
% Reorder dcm_infSA and ser_ind according to position_z
ind_in = ind_in(end:-1:1);
position_z = position_z(ind_in);
ser_ind = ser_ind(ind_in);
dcm_infSA = dcm_infSA(ind_in);
NZ_SA = length(ind_in);  % may have reduced number of SA slices

% Search for matching struct in setstruct
setstruct_ind = -1;
for i = 1:length(setstruct)
    if(sum(abs(orientation_1 - setstruct(i).ImageOrientation')) < IS_ZERO)
        setstruct_ind = i;
    end
end

if(setstruct_ind == -1)
    error('ConvertMaskToStack: error setting setstruct_ind Could not match segment file to image')
end

% Find first slice that matches 
first_segment_slice  = -1;
for i = 1:NZ_SA
    if(sum(abs(setstruct(setstruct_ind).ImagePosition' - dcm_infSA{i}.ImagePositionPatient)) < IS_ZERO)
        first_segment_slice = i;
        break;
    end
end
if(first_segment_slice == -1)
    error('Could not find slice matching first segment slice')
end
% NZ_Segment = setstruct(setstruct_ind).ZSize;

%% Write the first image from each series to a new series UID with instance
SeriesNumber = 1002;
SeriesUID = dicomuid;
SeriesDescription = 'SA Stack 3D';
for i = 1:NZ_SA
    dcm_im    = dicomread(dcm_infSA{i}.Filename);
    
    dcm_infSA{i}.SeriesNumber = SeriesNumber;
    dcm_infSA{i}.SeriesInstanceUID = SeriesUID;
    dcm_infSA{i}.SOPInstanceUID = dicomuid;
    dcm_infSA{i}.SeriesDescription = SeriesDescription;
    dcm_infSA{i}.MRAcquisitionType = '3D';
    dcm_infSA{i}.CardiacNumberOfImages = 1;
    
    fname_out = [dir_out,'write_',dcm_infSA{i}.SOPInstanceUID,'.dcm'];
    dicomwrite(dcm_im,fname_out,dcm_infSA{i},'WritePrivate',true,'CreateMode','Copy');
end
    
%% Generate 3D mask
% determin which setstruct to use based on series discription
% setstruct_ind set above
NX = setstruct(setstruct_ind).XSize;
NY = setstruct(setstruct_ind).YSize;
NZ = setstruct(setstruct_ind).ZSize;

sf = 7;
mask3D = zeros(NX*sf,NY*sf,NZ);

for i = 1:NZ
    mask = zeros(NX*sf,NY*sf);
    if ~isnan(setstruct(setstruct_ind).EndoX(1,1,i))
        mask = double(roipoly(mask,setstruct(setstruct_ind).EndoY(:,1,i)*sf,setstruct(setstruct_ind).EndoX(:,1,i)*sf));
    end
    mask3D(:,:,i) = mask;
end

% blockPlot(mask3D);

%% Correct for breathhold to breatthold variations
% Run xcorr2D to find slice to slice relationships, 
% then find slice most likely at end experation, 
% shift slices relative to that one, choosing smalles shift options.
breath_dir_z = [orientation_1(3) orientation_1(6)];
mag_bdz = sqrt(sum(breath_dir_z.^2));
shifts_xy = zeros(NZ,2);
corr_xy = zeros(NZ,2);
for i = 2:NZ
    if(sum(sum(mask3D(:,:,i-1))) ~= 0) && (sum(sum(mask3D(:,:,i))) ~= 0)
        tt = xcorr2(mask3D(:,:,i),mask3D(:,:,i-1));
        [tx ty] = ind2sub(size(tt),find((tt==max(max(tt)))));
        tx = tx -NX*sf; ty = ty-NY*sf;
        shift_ind = ([ty, tx]./repmat(sqrt(sum([ty, tx].^2,2)),1,2))*breath_dir_z'; % measure of how this corrleatets with Z
        [a ind] = max(abs(shift_ind));
        % Take component that agrees with breath_dir_z only
        mag_shift = sqrt(sum(([ty(ind), tx(ind)]/sf).^2));
        shifts_xy(i,:) = mag_shift*shift_ind(ind)*breath_dir_z/mag_bdz;  % only use component aligned with Z
        corr_xy(i,:) = [ty(ind), tx(ind)]/sf;
    end
end
%use the z component of orientation_1(3) and (6) to creat a vector of restricted direction
%subtract from tx and ty NX and NY to determin offset

%% find slice at furtherst exhale
shift_xy_cum = cumsum(shifts_xy);
if(reference_slice == -1)
    shift_ind = shift_xy_cum*breath_dir_z'; 
    [a ind] = max(shift_ind);
else
    ind = reference_slice;
end
shift_xy_cum = shift_xy_cum - repmat(shift_xy_cum(ind,:),NZ,1);

%%
setstruct_mod = setstruct;
[R1 R2 R3] = size(setstruct(setstruct_ind).EndoX);
setstruct_mod(setstruct_ind).EndoX = setstruct(setstruct_ind).EndoX - permute(repmat(shift_xy_cum(:,2),[1 R2,R1]),[3 2 1]);
setstruct_mod(setstruct_ind).EndoY = setstruct(setstruct_ind).EndoY - permute(repmat(shift_xy_cum(:,1),[1 R2,R1]),[3 2 1]);

%%
% Save the struct back to the file.
setstruct_orig = setstruct;
setstruct = setstruct_mod;
save(segment_out,'im', 'info', 'preview', 'setstruct');
setstruct = setstruct_orig;
clear('setstruct_orig');
%%
mask3D_im = zeros(NX,NY,NZ);

for i = 1:NZ
    mask = zeros(NX,NY);
    if ~isnan(setstruct_mod(setstruct_ind).EndoX(1,1,i))
        mask = double(roipoly(mask,setstruct_mod(setstruct_ind).EndoY(:,1,i),setstruct_mod(setstruct_ind).EndoX(:,1,i)));
    end
    mask3D_im(:,:,i) = mask;
end
figure;
blockPlot(mask3D_im);
title('Mask after alignment')

%% Write the segment mask from each series to a new series UID with instance
%load('G:\Data\4DFlow_WholeHeart\DCMIV005\DCMIV005_valign.mat')
SeriesNumber = 1003;
SeriesUID = dicomuid;
SeriesDescription = 'Segment Mask';

for i = 1:length(ser_ind)
    mask = uint16(mask3D_im(:,:,i)*100);
    
    dcm_infSA{i}.SeriesNumber = SeriesNumber;
    dcm_infSA{i}.SeriesInstanceUID = SeriesUID;
    dcm_infSA{i}.SOPInstanceUID = dicomuid;
    dcm_infSA{i}.SeriesDescription = SeriesDescription;
    dcm_infSA{i}.MRAcquisitionType = '3D';
    dcm_infSA{i}.CardiacNumberOfImages = 1;
    dcm_infSA{i}.SequenceName = 'Mask';
    
    fname_out = [dir_out, filesep,'write_',dcm_infSA{i}.SOPInstanceUID,'.dcm'];
    dicomwrite(mask,fname_out,dcm_infSA{i},'WritePrivate',true,'CreateMode','Copy');
end