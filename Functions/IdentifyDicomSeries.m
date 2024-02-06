function app = IdentifyDicomSeries(app,dicomData)
% Based on the function by Aaron Hess

% JK removed to compile app
if ~isdeployed
    [lib_dir,~,~]= fileparts(mfilename('fullpath'));
    cd(lib_dir);
    addpath(genpath([lib_dir,filesep,'Functions']));
end

NSer = length(dicomData.study(1).series);
name_list = cell(NSer,1);
ser_list = zeros(NSer,1);
inst_list = zeros(NSer,1);
is_phase = zeros(NSer,1);

for iSer = 1:NSer
    s = dicomData.study(1).series(iSer);
    name_list{iSer} = s.SeriesDescription;
    ser_list(iSer) = s.SeriesNumber;
    inst_list(iSer) = length(s.instance);
    
    try
        info = dicominfo(s.instance(1).Filename);
    catch
        uialert(app.OCMR4DFlowPostProcessingToolUIFigure,'Could not read DICOM information.','An error occured');
        app.DialogBoxTextArea.Value{length(app.DialogBoxTextArea.Value)+1} = char(string(datetime('now','Format','HH:mm')) + ' - Could not read DICOM information.'); pause(0.1); scroll(app.DialogBoxTextArea, 'bottom');
        error('Could not read DICOM information.');
    end
    is_phase(iSer) = isfield(info,'ImageType')&&((~isempty(strfind(info.ImageType,'\P\'))) || strcmp(info.ImageType(end-1:end),'\P'));
end

[ser_list,ind] = sort(ser_list);
valueIdx=1;
valuesStr={};
for iSer = 1:NSer
    % Update GUI
    if(is_phase(ind(iSer)))
        app.ListBox.Items{iSer} = [num2str(ser_list(iSer)),' : ',name_list{ind(iSer)},' : Phase'];
    else
        app.ListBox.Items{iSer} = [num2str(ser_list(iSer)),' : ',name_list{ind(iSer)},' : Magnitude'];
    end
    app.ListBox.ItemsData{iSer} = app.ListBox.Items{iSer};
    
    % Decide to highlight pre-selected flow series IDs
    if (any(ser_list(iSer) == app.serFlow))
        if(is_phase(ind(iSer)))
            valuesStr{valueIdx} = [num2str(ser_list(iSer)),' : ',name_list{ind(iSer)},' : Phase'];
        else
            valuesStr{valueIdx} = [num2str(ser_list(iSer)),' : ',name_list{ind(iSer)},' : Magnitude'];
        end
        valueIdx=valueIdx+1;
    end
end

% Show selected flow series
app.ListBox.Value={};
app.ListBox.Multiselect="on";
app.ListBox.Value=valuesStr;
end
