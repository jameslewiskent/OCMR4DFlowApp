function [app,ser_flow,dicomData] = GetFlowSeries(app)
% Get list of flow series from dicom folder.
% The first index is a magnitude series and remaining three are phase series (_P)
% All must have the same number of instances.
% Currently, this will only return the flow series for the scan on which the contouring was performed.
% If more than one flow scan was performed in a study, then this will be
% ignored, need to update code to fix this.
% JK Feb 2024
d_path = app.directoryPath;
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Extracting contours',"Indeterminate",'on'); pause(0.1)
app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Extracting CMR42 contours.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
try
    [contours, ~] = ExtractCMR42Contours(fullfile(d_path,app.ContourDropDown.Value)); % Find magnitude used for contouring
catch
    uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'An error occured extracting contours. ','An error occured','icon','error');
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - An error occured! '); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    error('An error occured extracting contours.')
end
close(progbar); % Below has its own progress bar

DicomFolderInfo = dir(d_path); DicomFolderNames = {DicomFolderInfo.name};
d_path = [d_path,filesep,DicomFolderNames{contains(DicomFolderNames,app.DefaultDicomDirName,'IgnoreCase',1)}];

ser_flow = zeros(1,4); % pre-allocate array

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

NStudies = size(dicomData.study,2); % Number of studies
for study_index = 1:NStudies
    NSeries(1,study_index) = size(dicomData.study(study_index).series,2); % Number of series in each study
end

series_total = 1; series_index = 1; study_index = 1; Mag_flag = 0; Study_Number = 1; clear ind is_match SOPUIDs
% First find the magnitude series ID 
 while series_total <= sum(NSeries) && Mag_flag ~= 1
    % Test contours were made on this magnitude image
    for ninstance = 1:size(dicomData.study(study_index).series(series_index).instance,2)
        SOPUIDs{ninstance} = dicomData.study(study_index).series(series_index).instance(ninstance).SOPInstanceUID;
    end
    for iImg = 1:length(contours)
         is_match = strfind(SOPUIDs,contours(iImg).uid);
         ind = find(~cellfun(@isempty,is_match), 1);
    end
    if ~isempty(ind)
        ser_flow(1,1) = dicomData.study(study_index).series(series_index).SeriesNumber;
        Mag_flag = 1; % We found the magnitude series number, yay!
        Study_Number = study_index;
    end
    
    % Increment counters
    series_index = series_index + 1;
    series_total = series_total + 1; % not same as series index in case of multiple studies
    if series_index > NSeries(1,study_index)
        study_index = study_index + 1; % increment study counter
        series_index = 1; % reset series index counter
    end
end

% Now we have found the magnitude dataset used for contouring, we will use
% the same study (Study_Number) for the remaining phase datasets
if Mag_flag == 0
    % Uh-Oh we still haven't found our magnitude dataset...
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Warning! Could not automatically find the magnitude series ID.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    Potential_Phase_Series_Indices= [];
else
    % Calculate approximate sequence duration
    Dicominfo = dicominfo(dicomData.study(Study_Number).series(series_index-1).instance(1).Filename);
    % JK Calculate approximate sequence timing from
    % dicom header, AFAiK there isn't a better way
    % of doing this since the actual sequence
    % time isn't stored in the header, so this is
    % approximate not accounting for recon time
    str1 = char(Dicominfo.SeriesTime); str2 = char(Dicominfo.AcquisitionTime);
    dhr = (str2double(str1(1:2))-str2double(str2(1:2)));
    dmin = (str2double(str1(3:4))-str2double(str2(3:4)));
    dsec = (str2double(str1(5:6))-str2double(str2(5:6)));
    mins = floor((dhr*3600 + dmin*60 + dsec)/60);
    secs = (dhr*3600 + dmin*60 + dsec) - mins * 60;
    app.total_seqtimemmss = [num2str(mins),':',num2str(secs)];
    app.total_seqtimeseconds = num2str(mins*60 + secs); 
    
    Mag_Img_Sz = size(dicomread(dicomData.study(Study_Number).series(find([dicomData.study(Study_Number).series.SeriesNumber] == ser_flow(1,1))).instance(1).Filename));
    
    % Now find all potential corresponding phase series IDs
    Phase_counter = 1; series_index = 1; study_index = 1;  Potential_Phase_Series_Indices = 0;
    for series_total = 1:sum(NSeries)
        
        Dicominfo = dicominfo(dicomData.study(study_index).series(series_index).instance(1).Filename);
        
        % First check whether SequenceName field exists in dicom, otherwise will throw error
        if isfield(Dicominfo,'SequenceName')
            % Phase series images must be the same size as magnitude images
            Current_Img_Sz = size(dicomread(dicomData.study(study_index).series(series_index).instance(1).Filename));
            if study_index == Study_Number && all(Current_Img_Sz == Mag_Img_Sz) && contains(Dicominfo.SeriesDescription,{'_P'}) && isfield(Dicominfo,'ImageType')&&((~isempty(strfind(Dicominfo.ImageType,'\P\'))) || strcmp(Dicominfo.ImageType(end-1:end),'\P'))
                Potential_Phase_Series_Indices(1,Phase_counter) = series_index; Phase_counter = Phase_counter + 1;
            end
        end
        
        series_index = series_index + 1;
        if series_index > NSeries(1,study_index)
            study_index = study_index + 1;     % increment study counter
            series_index = 1; % reset series index counter
        end
    end
end

% If there are only exactly 3 phase series indices, job done!
if size(Potential_Phase_Series_Indices,2) == 3
    for Phase_n = 1:3
        ser_flow(1,Phase_n+1) = dicomData.study(Study_Number).series(Potential_Phase_Series_Indices(1,Phase_n)).SeriesNumber;
    end
    
elseif size(Potential_Phase_Series_Indices,2) > 3
    % Else we must choose which phase series indices we use.
    % Currently, do this based on closest series time to the
    % magnitude
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - WARNING: Multiple potential phase series flow IDs found. Choosing based on nearest in time to magntiude flow series ID. Please double check the chosen series IDs.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    disp(['Warning: @',mfilename,': Multiple potential phase series flow IDs found. Choosing based on nearest in time to magntiude flow series ID. Please double check the chosen series IDs.'])
    app.SeriesIDsEditField.FontColor = [1 0.5 0];
    
    Mag_SeriesTime = str2double(dicomData.study(Study_Number).series(find([dicomData.study(Study_Number).series.SeriesNumber] == ser_flow(1,1))).SeriesTime);
    for Phase_n = 1:size(Potential_Phase_Series_Indices,2)
        Time_Diff(1,Phase_n) = abs(str2double(dicomData.study(Study_Number).series(Potential_Phase_Series_Indices(1,Phase_n)).SeriesTime) - Mag_SeriesTime);
    end
    % Find three smallest indices
    [~,min_index] = mink(Time_Diff,3,2);
    Potential_Phase_Series_Indices = Potential_Phase_Series_Indices(1,sort(min_index));
    for Phase_n = 1:3
        ser_flow(1,Phase_n+1) = dicomData.study(Study_Number).series(Potential_Phase_Series_Indices(1,Phase_n)).SeriesNumber;
    end
end

% Check we found all the series indices, if any indices are still 0, we
% have a problem!
if ser_flow(1,1) == 0
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - WARNING: Failed to find a dicom series which corresponds to the supplied contour file!'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    disp(['Warning: @',mfilename,': Failed to find a dicom series that corresponds to the supplied contour file.'])
    app.SeriesIDsEditField.FontColor = [1 0 0];
elseif any(ser_flow == 0)
    disp(['Warning: @',mfilename,': Failed to find all dicom flow series! Please enter series IDs manually in preferences tab.'])
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - WARNING: Failed to find all dicom flow series! Please enter series IDs manually in preferences tab.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    app.SeriesIDsEditField.FontColor = [1 0 0];
elseif any(ser_flow(1) == ser_flow(2:end))
    % In this scenario, we are going to reassign the first
    % index to a magnitude dicom flow series which
    % should be the next flow series ID
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Expected contouring was performed on the magnitude image, but it appears to have been done on the phase image.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    app.CountoursonPhaseImage = ser_flow(1);
    ser_flow(1) = ser_flow(1) + 1; % Currently this is a phase series ID, next one should be a magnitude
end
app.serFlow = ser_flow;
app.StudyNumber = Study_Number;
app.DICOMStudyEditField.Value = Study_Number; % Show on UI
app.DICOMStudyEditField.Editable = 'off'; % Prevent editing is automatically found

close(progbar)
end

