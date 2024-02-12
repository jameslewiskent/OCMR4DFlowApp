function app = finddicoms(app)
% Load dicoms and find series IDs
    app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Begining to load in DICOM files and contours. Please be patient the next few steps may take a few minutes.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
    
    try       
        % Check Dicoms are in correct format and dicoms are in folders
        % named by their series name
        d = dir(app.directoryPath);
        for i = 1:numel(d)
            if contains(string(d(i).name),string(app.DefaultDicomDirName),'IgnoreCase',1)
                dicomfolderpresent = 1; break
            else
                dicomfolderpresent = 0;
            end
        end
        if ~dicomfolderpresent  % if no dicom folder reformat and put into DICOM folder
            app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Files not the in expected DICOM format. Converting to DICOM.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
            Changetodcm(app.directoryPath,app);
        end
        d = dir(fullfile(app.directoryPath,app.DefaultDicomDirName));
        d([d.isdir])= []; % Remove all directories.

        % Ignore .ima.mat files
        for n = 1:length(d)
           if contains(d(n).name,'.ima.mat','IgnoreCase',true)
               ignore_list(n) = true;
           else
               ignore_list(n) = false;
           end
        end
        if exist('ignore_list','var')
        d(ignore_list) = [];
        end
        
        if ~isempty(d) % if dicoms not in folders named by their series
            app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - DICOMs not in the expected file structure. Re-structuring.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
            DicomMoveToSeriesDir(fullfile(app.directoryPath,app.DefaultDicomDirName),app);
        end
        
        [~,app.CurrentDirectoryEditField.Value] = fileparts(app.directoryPath);
        app.CurrentDirectoryEditField.Tooltip = app.directoryPath;
        
        try
            [app,ser_flow,dicomData] = GetFlowSeries(app); pause(0.1) % JK select a default flow series based on dicom series description / sequence
        catch
            dicomData = processDicomDirRecursive(app,app.directoryPath,'*'); % If above has failed dicomData doesn't exist
            ser_flow = app.SeriesIDsEditField.Value; % Use default values
            app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Could not automatically find flow series IDs. Please enter these manually.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
            error('Could not automatically find flow series IDs.')
        end
      
        Series_String = sprintf('%.0f,' , ser_flow);
        Series_String = ['[',Series_String(1:end-1),']']; % strip final comma, encase in square brackets
        app.SeriesIDsEditField.Value = Series_String; % JK update string in preferences with the selected flow series
        
        app.ListBox.Items = {}; % JK clear the list of items in listbox when new path specified
        
        if length(dicomData.study) > 1
        app.DICOMStudyEditField.Limits = [1,length(dicomData.study)]; % Limit to number of studies
        else
        app.DICOMStudyEditField.Editable = 'off';    
        end
        
        IdentifyDicomSeries(app,dicomData); % Update selection
        
        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Succesfully loaded dicoms. Ready for analysis.'); scroll(app.DialogBoxTextArea, 'bottom');pause(0.1)
    catch % else if failed
        uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'An error occured and the dicoms could not be loaded.','An error occured','icon','error');
        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - An error occured and the dicoms could not be loaded.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
        app.GoButton.Text = 'Error Occured!';
        app.Lamp.Color = [1 0 0];
    end
end

