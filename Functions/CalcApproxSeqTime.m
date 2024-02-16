function [mins,secs] = CalcApproxSeqTime(dicomData,Study_Number,Series_Index)
% Calculate approximate sequence duration
% JK Calculate approximate sequence timing from
% dicom header, AFAIK there isn't a better way
% of doing this since the actual sequence
% time isn't stored in the header, so this is
% approximate not accounting for recon time

Dicominfo = dicominfo(dicomData.study(Study_Number).series(Series_Index).instance(1).Filename);

str1 = char(Dicominfo.SeriesTime); str2 = char(Dicominfo.AcquisitionTime);
dhr = (str2double(str1(1:2))-str2double(str2(1:2)));
dmin = (str2double(str1(3:4))-str2double(str2(3:4)));
dsec = (str2double(str1(5:6))-str2double(str2(5:6)));
mins = floor((dhr*3600 + dmin*60 + dsec)/60);
secs = (dhr*3600 + dmin*60 + dsec) - mins * 60;
end

