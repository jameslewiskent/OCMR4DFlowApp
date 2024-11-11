function [app,ser_flow,dicomData] = GetFlowSeries(app)
% Get the list of flow series IDs from dicom folder.
% All must have the same number of instances.
% Currently, this will only return the flow series for the scan on which the contouring was performed.
% JK updated Feb 2024

% Update UI progress bar
d_path = app.directoryPath;
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Extracting contours',"Indeterminate",'on'); pause(0.1)
app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Extracting CMR42 contours.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');

% Try to get the study and image UID which were used for contouring
try
    [~, study_uid, ~, image_uid] = ExtractCMR42Contours(fullfile(d_path,app.ContourDropDown.Value)); % Find magnitude used for contouring
catch
    uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'An error occured extracting contours. ','An error occured','icon','error');
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - An error occured! '); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    error('An error occured extracting contours.')
end
close(progbar); % Below has its own progress bar

DicomFolderInfo = dir(d_path); DicomFolderNames = {DicomFolderInfo.name};
d_path = [d_path,filesep,DicomFolderNames{contains(DicomFolderNames,app.DefaultDicomDirName,'IgnoreCase',1)}];

ser_flow = zeros(1,4); % pre-allocate array

% Read in DICOM data and update progress bar
app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Reading in DICOM files.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
try
    dicomData = processDicomDirRecursive(app,d_path,'*');
catch
    uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'An error occured processing DICOMs recursive.','An error occured','icon','error');
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - An error occured!. '); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    error('An error occured processing DICOMs recursive.');
end
app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Attempting to detect flow series IDs.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Finding DICOM Flow Series IDs',"Indeterminate",'on'); pause(0.1)

% First, if there are multiple studies, let's find the study number the contours are associated with
if size(dicomData.study,2) >1
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Multiple studies found inside DICOM data. Selecting study corresponding to contours. '); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    for study_index = 1:size(dicomData.study,2)
        if strcmp(dicomData.study(study_index).StudyInstanceUID,study_uid)
            Study_Number = study_index;
        end
    end
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + [' - DICOM Study: ',num2str(Study_Number)]); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
else
    Study_Number = 1;
end

NSeries = size(dicomData.study(Study_Number).series,2); % Number of series in the chosen study

% Now let's find the series ID which is associated with the contouring
series_index = 1; flag = 0;
while series_index <= NSeries && flag ~= 1
    for ninstance = 1:size(dicomData.study(Study_Number).series(series_index).instance,2)
        SOPUIDs{ninstance} = dicomData.study(Study_Number).series(series_index).instance(ninstance).SOPInstanceUID;
    end
    if any(contains(SOPUIDs,image_uid))
        Contour_Series_ID_Found = dicomData.study(Study_Number).series(series_index).SeriesNumber;
        Series_Indices(1,1) = series_index;
        ser_flow(1,1) = Contour_Series_ID_Found;
        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + [' - Contouring performed on series ID: ',num2str(ser_flow(1,1))]); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
        flag = 1; % We found a series number, yay!
    end
    % Increment counter
    series_index = series_index + 1;
end

if flag == 0
    % Uh-Oh we still haven't found our magnitude dataset...
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - WARNING: Failed to find all dicom flow series! Please enter series IDs manually in preferences tab.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
end

% Find ProtocolName which the series ID we found belows to
seriesID_ProtocolName = dicominfo(dicomData.study(Study_Number).series(Series_Indices(1,1)).instance(1).Filename).ProtocolName;

% Find all series which have the same protocol name
n = 2;
for series_index = 1:NSeries
    if isfield(dicominfo(dicomData.study(Study_Number).series(series_index).instance(1).Filename),'ProtocolName') % Check field exists
        if strcmp(seriesID_ProtocolName,dicominfo(dicomData.study(Study_Number).series(series_index).instance(1).Filename).ProtocolName)
            ser_flow(1,n) = dicominfo(dicomData.study(Study_Number).series(series_index).instance(1).Filename).SeriesNumber;
            Series_Indices(1,n) = series_index;
            n = n + 1;
        end
    end
end

% remove duplicates but keep order! (DO NOT USE UNIQUE)
[bs, vec] = sort(ser_flow(:).'); uvec(vec) = [true, diff(bs) ~= 0]; ser_flow   = ser_flow(uvec);
[bs, vec] = sort(Series_Indices(:).'); uvec(vec) = [true, diff(bs) ~= 0]; Series_Indices = Series_Indices(uvec);

% Now make sure all have same dimensions and find image type
for n = 1:length(Series_Indices)
    Img_Sz(:,n) = size(dicomData.study(Study_Number).series(Series_Indices(n)).instance,2);
    info = dicominfo(dicomData.study(Study_Number).series(Series_Indices(n)).instance(1).Filename);
    Is_Phase_Img(n) = ((~isempty(strfind(info.ImageType,'\P\'))) || strcmp(info.ImageType(end-1:end),'\P'));
end

% Remove those without most common size
for n = 1:length(Series_Indices)
    Remove_flag(:,n) = size(dicomData.study(Study_Number).series(Series_Indices(n)).instance,2) ~= mode(Img_Sz);
end
Series_Indices(Remove_flag) = [];
ser_flow(Remove_flag) = [];
Is_Phase_Img(Remove_flag) = [];

% If the contour series ID we originally found is a subset of the
% remaining group and a magnitude image, use this and the three phases found
if ~Is_Phase_Img(ser_flow == Contour_Series_ID_Found)
    
    if length(ser_flow(Is_Phase_Img)) > 3
        % There are more than 3 phases to choose from still!
        % For now, let's choose those with the closest series time to the magnitude series
        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - WARNING: Multiple potential phase series flow IDs found. Choosing based on nearest in time to magntiude flow series ID. Please double check the chosen series IDs.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
        disp(['Warning: @',mfilename,': Multiple potential phase series flow IDs found. Choosing based on nearest in time to magntiude flow series ID. Please double check the chosen series IDs.'])
        app.SeriesIDsEditField.FontColor = [1 0.5 0];
        
        Potential_Phase_Series_Indices = Series_Indices(Is_Phase_Img);
        
        Mag_SeriesTime = str2double(dicomData.study(Study_Number).series(find([dicomData.study(Study_Number).series.SeriesNumber] == Contour_Series_ID_Found)).SeriesTime);
        for Phase_n = 1:size(Potential_Phase_Series_Indices,2)
            Time_Diff(1,Phase_n) = abs(str2double(dicomData.study(Study_Number).series(Potential_Phase_Series_Indices(1,Phase_n)).SeriesTime) - Mag_SeriesTime);
        end
        
        % Find three smallest indices
        [~,min_index] = mink(Time_Diff,3,2);
        Potential_Phase_Series_Indices = Potential_Phase_Series_Indices(1,sort(min_index));
        for Phase_n = 1:3
            ser_flow(1,Phase_n+1) = dicomData.study(Study_Number).series(Potential_Phase_Series_Indices(1,Phase_n)).SeriesNumber;
        end
        
    else
        ser_flow(1,2:4) = ser_flow(Is_Phase_Img);
    end
    ser_flow(1,1) = Contour_Series_ID_Found;
    ser_flow = ser_flow(1,1:4); % clear remaining from selection
end

% Calculate the appproximate sequence time based on dicom information
for n = 1:length(Series_Indices)
    [mins(n),secs(n)] = CalcApproxSeqTime(dicomData,Study_Number,Series_Indices(n));
end
mins = round(mean(mins)); secs = round(mean(secs));
app.total_seqtimemmss = [num2str(mins),':',num2str(secs)];
app.total_seqtimeseconds = num2str(mins*60 + secs);

% Check we found all the series indices, if any indices are still 0, we have a problem!
if any(ser_flow == 0)
    disp(['Warning: @',mfilename,': Failed to find all dicom flow series! Please enter series IDs manually in preferences tab.'])
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - WARNING: Failed to find all dicom flow series! Please enter series IDs manually in preferences tab.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    app.SeriesIDsEditField.FontColor = [1 0 0];
end
app.serFlow = ser_flow;
app.StudyNumber = Study_Number;
app.DICOMStudyEditField.Value = Study_Number; % Show on UI
app.DICOMStudyEditField.Editable = 'off'; % Prevent editing is automatically found

close(progbar)
end

