% Enter the directory to search
directory = 'C:\Users\antonios\Documents\TAVI\TAVI_05\DICOM\12420004';
% List all .out files
fileList = dir([directory, '\*']); 
% Loop through each .out file, copy it and give new extension: .txt
for i = 1:numel(fileList)
    file = fullfile(directory, fileList(i).name);
    [tempDir, tempFile] = fileparts(file); 
    status = copyfile(file, fullfile(tempDir, [tempFile, '.dcm']));
    delete(file);
end