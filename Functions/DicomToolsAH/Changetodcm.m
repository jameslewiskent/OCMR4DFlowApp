function Changetodcm(d_path,app)
% List all .out files
fileList = dir(fullfile(d_path, filesep,'*')); 
fileList([fileList.isdir])= []; % JK temp remove directories... This was causing a recursive issue
for i = 1:numel(fileList)
    if contains(string(fileList(i).name),string('cvi42ws'),'IgnoreCase',1)
fileList(i)= []; % JK temp remove any contour file potentially in the same folder
    end
end

mkdir(fullfile(d_path,'DICOM')); % JK added make new DICOM directory

% Loop through each .out file, copy it and give new extension: .dcm
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please wait','Message','Changing input format to DICOM');
for i = 1:numel(fileList)
    file = fullfile(d_path, fileList(i).name);
    [tempDir, tempFile] = fileparts(file); 
    copyfile(file, fullfile([tempDir,filesep,'DICOM'], [tempFile, '.dcm'])); % Copy to the new DICOM directory
    delete(file);
  
    progbar.Value = i/numel(fileList);
end
progbar.Value = 1;
close(progbar);
end

