function dicomData = processDicomDirRecursive(app,strPath, varargin)
% Return dicomData structure for all DICOM files in a folder recursively.
%
% strPath - folder to search.
%
% Any further arguments (e.g. wildcard) are passed unchanged to
% processDicomDir().
%
% Example:
% dicomData = processDicomDirRecursive(pwd)

% Copyright Chris Rodgers, University of Oxford, 2011.
% $Id$

% Perform a depth-first traversal of folders under strPath
dirList{1} = strPath;

idx = 1;
while idx <= numel(dirList)
    d = dir(dirList{idx});
    
    d(strcmp({d.name},'.')) = [];
    d(strcmp({d.name},'..')) = [];
    
    newDirs = {d([d.isdir]).name}.';
    
    for jdx=1:numel(newDirs)
       newDirs{jdx} = fullfile(dirList{idx},newDirs{jdx}); 
    end
    
    dirList = [dirList(1:idx); newDirs; dirList((idx+1):end)];
    
    idx = idx+1;
end

% JK moved progress bar and incorporated into UI figure window
progbar = uiprogressdlg(app.OCMR4DFlowPostProcessingToolUIFigure,'Title','Please Wait','Message','Loading DICOM files'); pause(0.1)
% Process each individual folder
for idx=1:numel(dirList)
    dicom{idx} = processDicomDir(dirList{idx},varargin{:});
    progbar.Value = idx./numel(dirList);
    progbar.Message = ['Loading DICOM file ',num2str(idx),' of ', num2str(numel(dirList))];
end
close(progbar);

% Merge the results into a single dicomData struct
dicomData = dicomDataMerge(dicom{:});
