function [dicomData] = processDicomDir(strDir, strWildcard)
% Scan all files in a directory, find the DICOM files and sort them.
%
% Optionally, restrict the search to files matching a wildcard.

% TODO: Make sorting this tree optional.
% TODO: Make progress bar optional.
% TODO: Allow choice of fields to be extracted.
% TODO: Deal with hash collisions in the cache.
% TODO: Deal with cache files that have different fields.
% TODO: Allow override of the cache filenames.

% $Id: processDicomDir.m 4592 2011-10-13 13:11:02Z crodgers $
% Copyright Chris Rodgers, University of Oxford, 2010-11.

% File format counter: increment to force refresh of cache files.
version = 4;

if nargin < 2
    strWildcard = '';
end

%% Check folder exists
if ~exist(strDir,'dir')
    error('RodgersSpectroTools:DirectoryNotFound','Directory not found ("%s").',strDir);
end

%% Check for cache file
cacheFile = fullfile(tempdir(),['processDicomDir_CACHE_' dirHash(strDir, strWildcard) '.mat']);
try
    cacheData = load(cacheFile);
    
    if cacheData.version == version
        dicomData = cacheData.dicomData;
        return
    else
        clear cacheData
    end
catch ME
    if ~strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
        rethrow(ME);
    end
end

%% No matching cache file was found, so we must actually scan directory.

% Initialise dcm4che2 toolkit
dcm4che2_init

% Scan for DICOM files
if nargin >= 2 && ~isempty(strWildcard)
    % Append wildcard if it was supplied
    strDirWithWildcard = fullfile(strDir,strWildcard);
    files = dir(strDirWithWildcard);
elseif exist(fullfile(strDir, 'DICOMDIR'),'file') == 2
    % This directory is a DICOMDIR tree. E.g. CD/DVD.
    % For now, extract all the filenames from the DICOMDIR and then analyse
    % each one in turn below.
    %
    % TODO: For speed, extract the desired parameters directly from the
    %       DICOMDIR.
    
    theDicomdir = dicominfo(fullfile(strDir, 'DICOMDIR'));
    
    fn = fieldnames(theDicomdir.DirectoryRecordSequence);
    files = [];
    for idx=1:numel(fn)
        if isfield(theDicomdir.DirectoryRecordSequence.(fn{idx}),'ReferencedFileID')
            % The regexp fixes path separators in the DICOMDIR
            files(end+1).name = regexprep(theDicomdir.DirectoryRecordSequence.(fn{idx}).ReferencedFileID, '[\\/]', regexptranslate('escape',filesep()));
            files(end).isdir = 0;
        end
    end
else
    files = dir(strDir);
end

dicomData.path = strDir;
dicomData.wildcard = strWildcard;
dicomData.study = struct('StudyInstanceUID',{},'StudyID',{},'StudyDescription',{},'series',{});

for idx = 1:numel(files)    
    thisFile = fullfile(strDir,files(idx).name);
    
    if ~files(idx).isdir
        % Process this DICOM file
        
        % Catch a Java error if it is not a DICOM file.
        % Profiling confirms that the try..catch block is faster than
        % calling dicomCheckMagic(...)
        try
            dcm = dcm4che2_readDicomFile(thisFile);
        catch ME
            % Catch errors that are due to this not being a DICOM file, but
            % report other errors e.g. out of memory.
            %
            % org.dcm4che2.io.DicomCodingException: Not a DICOM Stream
            
            if strcmp(ME.identifier,'MATLAB:Java:GenericException')
                msg = regexp(ME.message,'^[^\n\r]*[\n\r]([^\n\r]*)','once','tokens');
                
                if strcmp(msg{1},'org.dcm4che2.io.DicomCodingException: Not a DICOM Stream') ...
                   || strcmp(msg{1},'java.io.EOFException')
                    continue;
                end
                
                % Ignore errors from non-DICOM files.
                % This check is slower for most files than that above.
                fprintf('\nError in file %s\n',thisFile);
                if ~dicomCheckMagic(thisFile)
                    continue;
                end
            end
           
            % Any other error should be rethrown
            fprintf('Error while processing file "%s":\n',thisFile)
            rethrow(ME);
        end
        
        % Build index of studies
        myStudyInstanceUID = char(dcm.getString(org.dcm4che2.data.Tag.StudyInstanceUID));
        
        myStudyDx = find(strcmp({dicomData.study.StudyInstanceUID},myStudyInstanceUID));
        if isempty(myStudyDx)
            myStudyDx = numel(dicomData.study) + 1;
            dicomData.study(myStudyDx).series = struct('SeriesInstanceUID',{},'SeriesNumber',{},'SeriesDescription',{},'instance',{});
        end
        dicomData.study(myStudyDx).StudyInstanceUID = myStudyInstanceUID;
        dicomData.study(myStudyDx).StudyID = char(dcm.getString(org.dcm4che2.data.Tag.StudyID));
        dicomData.study(myStudyDx).StudyDescription = char(dcm.getString(org.dcm4che2.data.Tag.StudyDescription));
        
        % Build index of series
        mySeriesInstanceUID = char(dcm.getString(org.dcm4che2.data.Tag.SeriesInstanceUID));
        
        mySeriesDx = find(strcmp({dicomData.study(myStudyDx).series.SeriesInstanceUID},mySeriesInstanceUID));
        if isempty(mySeriesDx)
            mySeriesDx = numel(dicomData.study(myStudyDx).series) + 1;
            dicomData.study(myStudyDx).series(mySeriesDx).instance = struct('SOPInstanceUID',{},'InstanceNumber',{});
        end
        dicomData.study(myStudyDx).series(mySeriesDx).SeriesInstanceUID = mySeriesInstanceUID;
        dicomData.study(myStudyDx).series(mySeriesDx).SeriesNumber = dcm.getInt(org.dcm4che2.data.Tag.SeriesNumber);
        dicomData.study(myStudyDx).series(mySeriesDx).SeriesDescription = char(dcm.getString(org.dcm4che2.data.Tag.SeriesDescription));
        dicomData.study(myStudyDx).series(mySeriesDx).SeriesDate = char(dcm.getString(org.dcm4che2.data.Tag.SeriesDate));
        dicomData.study(myStudyDx).series(mySeriesDx).SeriesTime = char(dcm.getString(org.dcm4che2.data.Tag.SeriesTime));
        
        % Build index of instances
        mySOPInstanceUID = char(dcm.getString(org.dcm4che2.data.Tag.SOPInstanceUID));
        
        % Modify so that all files present are loaded, even duplicates.
%         myInstanceDx = find(strcmp({dicomData.study(myStudyDx).series(mySeriesDx).instance.SOPInstanceUID},mySOPInstanceUID));
%         if isempty(myInstanceDx)
        myInstanceDx = numel(dicomData.study(myStudyDx).series(mySeriesDx).instance) + 1;
%         end
        dicomData.study(myStudyDx).series(mySeriesDx).instance(myInstanceDx).SOPInstanceUID = mySOPInstanceUID;
        dicomData.study(myStudyDx).series(mySeriesDx).instance(myInstanceDx).InstanceNumber = dcm.getInt(org.dcm4che2.data.Tag.InstanceNumber);
        dicomData.study(myStudyDx).series(mySeriesDx).instance(myInstanceDx).Filename = thisFile;
        dicomData.study(myStudyDx).series(mySeriesDx).instance(myInstanceDx).ImageComment = char(dcm.getString(org.dcm4che2.data.Tag.ImageComments));
    end
end

% Sort the data before returning because we'll almost always want this
for StudyDx=1:numel(dicomData.study)
    for SeriesDx=1:numel(dicomData.study(StudyDx).series)
        [tmp, sortDx] = sort([dicomData.study(StudyDx).series(SeriesDx).instance.InstanceNumber]);
        dicomData.study(StudyDx).series(SeriesDx).instance = dicomData.study(StudyDx).series(SeriesDx).instance(sortDx);
    end
    [tmp, sortDx] = sort([dicomData.study(StudyDx).series.SeriesNumber]);
    dicomData.study(StudyDx).series = dicomData.study(StudyDx).series(sortDx);
end
[tmp, sortDx] = sort({dicomData.study.StudyID}); % These are strings - need different sorting.
dicomData.study = dicomData.study(sortDx);

% Save in cache file.
try
    save(cacheFile,'dicomData','version');
catch
    warning('Cannot save cache file: "%s".', cacheFile);
end

end
        
function [dcm] = dcm4che2_readDicomFile(strFile)
% Load DICOM file with the dcm4che2 toolkit

% Uncomment these lines to profile DICOM file scanning.
% tic
din = org.dcm4che2.io.DicomInputStream(java.io.File(strFile));
dcm = din.readDicomObject();
% fprintf('%g ms for\t''%s''\n',toc(),strFile);

end

function dcm4che2_init()
% Check first that the Java path is set properly
% Then add the java libraries to your path
checkjava = which('org.dcm4che2.io.DicomInputStream');
if isempty(checkjava)
    libpath = fullfile(fileparts(mfilename('fullpath')),'dcm4che','dcm4che-2.0.24-bin','dcm4che-2.0.24','lib');
    %fprintf('\nlibpath = %s',libpath);
    javaaddpath(fullfile(libpath,'dcm4che-core-2.0.24.jar'));
    javaaddpath(fullfile(libpath,'dcm4che-image-2.0.24.jar'));
    javaaddpath(fullfile(libpath,'dcm4che-imageio-2.0.24.jar'));
    javaaddpath(fullfile(libpath,'dcm4che-iod-2.0.24.jar'));
    javaaddpath(fullfile(libpath,'slf4j-api-1.6.1.jar'));
    javaaddpath(fullfile(libpath,'slf4j-log4j12-1.6.1.jar'));
    javaaddpath(fullfile(libpath,'log4j-1.2.16.jar'));
    
    checkjava = which('org.dcm4che2.io.DicomInputStream');
    
    if isempty(checkjava)
        error('Cannot load dcm4che2 v2.0.24 toolkit.')
    end
end
end
