function PrintSeriesInfo(data_dir)
% Aaron Hess
% university of oxford
% Print dicom information including series numbers
    

    %[lib_dir,~,~]= fileparts(mfilename('fullpath'));
    %addpath([lib_dir,filesep,'Rodgers_dcm_read']);

    [dicomData] = processDicomDirRecursive(data_dir, '*'); 

    NStudies = length(dicomData.study);
    for iStudy = 1:NStudies
    NSer = length(dicomData.study(iStudy).series);
    name_list = cell(NSer,1);    
    ser_list = zeros(NSer,1);
    inst_list =zeros(NSer,1);
    is_phase = zeros(NSer,1);
    for iSer = 1:NSer
        s = dicomData.study(iStudy).series(iSer);
        name_list{iSer} = s.SeriesDescription;    
        ser_list(iSer) = s.SeriesNumber;
        inst_list(iSer) = length(s.instance);
        
        info = dicominfo(s.instance(1).Filename);
        is_phase(iSer) = isfield(info,'ImageType')&&((~isempty(strfind(info.ImageType,'\P\'))) || strcmp(info.ImageType(end-1:end),'\P'));
    end
    fprintf('\nDICOM SERIES INFORMATION\n');
    fprintf('SER   :           DESCRPTION               : IMAGE COUNT\n');
    [ser_list,ind] = sort(ser_list);
    for iSer = 1:NSer
        if(is_phase(ind(iSer)))
            fprintf('%5d : %35s : %5d : Phase\n',ser_list(iSer),name_list{ind(iSer)},inst_list(ind(iSer)));
        else
            fprintf('%5d : %35s : %5d : Magnitude\n',ser_list(iSer),name_list{ind(iSer)},inst_list(ind(iSer)));
        end
    end
    end
end