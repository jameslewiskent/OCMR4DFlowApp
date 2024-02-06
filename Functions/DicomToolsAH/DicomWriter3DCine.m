function DicomWriter3DCine(out_dir,im3D,is_bart,hdr,base_info,name)
% write a 3D cine to a new dicom series
% out_dir is where to place the files which will be named according to uid
% im3D is a 4D matrix ROxPE1xPE2xT 
% is_bart reflects the axes due to BART using ifft
% hdr is the hdr returend from mapvbvd with orientation infomration
% dcm_info is a starting dicom info from the same subject as determined by matlab

RO = 1;
PE1 = 2;
PE2 = 3;

TO_SLICE = RO;


% What fields must be filled in in dicom
info = base_info;  %start with one that works...

info.ColorType = 'grayscale';

info.StudyDescription = base_info.StudyDescription;

info.AcquisitionDate = base_info.AcquisitionDate;
info.AcquisitionTime = 0;
info.ContentDate = base_info.ContentDate;
info.ContentTime = base_info.ContentTime; % not sure how to find this
info.DeviceSerialNumber = num2str(hdr.Dicom.DeviceSerialNumber);
info.EchoNumber = 1;
info.EchoTime = hdr.Meas.alTE(1)/1000;
info.EchoTrainLength = 1;
info.FlipAngle = hdr.Meas.FlipAngle;
info.FrameOfReferenceUID = base_info.FrameOfReferenceUID;
info.ImagedNucleus = '1H';
info.ImagingFrequency = hdr.Config.SystemFrequency/1e6;
info.ImageType =  'USER';

info.InstitutionName = hdr.Dicom.InstitutionName;
info.InstitutionAddress = 'Headley Way StreetNo,Oxford/228790/,Oxfordshire,GB,OX3 9DU';

info.InstanceCreationDate =sprintf('%s', datetime('today','Format','yyyyMMdd'));
info.InstanceCreationTime = sprintf('%s',datetime(datetime,'Format','HHmmss'));

info.Modality = hdr.Dicom.Modality;
info.Manufacturer = hdr.Dicom.Manufacturer;
info.MagneticFieldStrength = hdr.Dicom.flMagneticFieldStrength;
info.ManufacturerModelName = 'Investigational_Device_7T';
info.NumberOfAverages = 1;
%info.PatientAge = 99;

info.PatientBirthDate= base_info.PatientBirthDate; % hdr.Config.PatientBirthDay;
info.PatientID = hdr.Config.PatientID;
info.PatientName =  hdr.Dicom.tPatientName;
info.PatientPosition = hdr.Dicom.tPatientPosition;
info.PatientSex = base_info.PatientSex;
info.PatientWeight = 50;
info.PercentPhaseFieldOfView = 100;
info.PercentSampling = 100; % not sure what this is
info.PixelBandwidth = 1e6/hdr.Meas.alDwellTime(1);
info.ProtocolName = hdr.Dicom.tProtocolName;
info.RepetitionTime = hdr.Meas.TR_us/1000;
info.SamplesPerPixel = 1;
info.SAR = 0;
info.ScanningSequence = hdr.Dicom.tScanningSequence;
info.SequenceName = hdr.Meas.SequenceString;
info.SeriesDate = sprintf('%s',datetime('today','Format','yyyyMMdd'));
info.SOPClassUID = base_info.SOPClassUID;
info.SpecificCharacterSet = base_info.SpecificCharacterSet;

info.SeriesDescription = [name,'_',hdr.Dicom.tProtocolName];

info.SeriesNumber = 1000 + round(10*rand(1));
info.SmallestImagePixelValue = 0;
% info.SoftwareVersion
info.StudyID = base_info.StudyID;
info.StudyDate = base_info.StudyDate;
info.StudyTime = base_info.StudyTime;
info.TransmitCoilName = hdr.Dicom.TransmittingCoil;

info.MRAcquisitionType = '3D';
info.LargestImagePixelValue = 4095;

% get from input example if there is one
info.SpacingBetweenSlices = 0;
info.SeriesInstanceUID = dicomuid;


% Work things out for this data
if(is_bart)
    im3D = flip(flip(flip(im3D,1),2),3);
end
% normalise range
im3D = round(abs(im3D)*4095)/max(abs(im3D(:)));

sz = size(im3D);

% resolution
meas.fov(1) = hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
meas.fov(2) = hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
meas.fov(3) = hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;

meas.matrix(1) = hdr.MeasYaps.sKSpace.lBaseResolution;
meas.matrix(2) = hdr.MeasYaps.sKSpace.lPhaseEncodingLines;
meas.matrix(3) = hdr.MeasYaps.sKSpace.lPartitions;

meas.pixel_size = meas.fov./meas.matrix;


% specify for this data

info.Height = sz(2);
info.Widt  = sz(3);
info.Rows = sz(2);
info.Columns = sz(3);

info.AcquisitionMatrix = sz(2:3);  % 2D

info.ImageOrientationPatient = [0 1 0 1 0 0];
info.ImagePositionPatient = [hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dSag, hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dCor, hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition.dTra];

info.ImagePositionPatient(3) = info.ImagePositionPatient(3) - (sz(1)/2)*meas.pixel_size(1);

info.InPlanePhaseEncodingDirection = 'ROW';
%info.
info.NumberOfPhaseEncodingSteps = sz(2);

info.PixelSpacing = meas.pixel_size(2:3);  % 2x1
info.SliceThickness = meas.pixel_size(1);
info.WindowCenter = 2048;
info.WindowWidth = 4095;

% for each time frame and slice
sl = info.ImagePositionPatient(3);
inst = 0;
for iSl = 1:sz(1)
    info.ImagePositionPatient(3) = sl + iSl*meas.pixel_size(1);
    info.SliceLocation = info.ImagePositionPatient(3);
    for iPh = 1:sz(4)
        inst = inst+1;
        info.InstanceNumber = inst;
        info.SOPInstanceUID = dicomuid;
        
        fname_out = [out_dir, filesep,'write_',num2str(iSl),'_',num2str(iPh),'_', info.SOPInstanceUID,'.dcm'];
        dicomwrite(uint16(squeeze(im3D(iSl,:,:,iPh))),fname_out,info,'WritePrivate',true,'CreateMode','Copy');
    end
end
end