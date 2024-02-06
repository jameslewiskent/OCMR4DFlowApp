function [contours, study_uid] = ExtractCMR42Contours(fname)
% Aaron Hess
% University of Oxford
% July 2015
% contact: aaron.hess@cardiov.ox.ac.uk
% If you use this function, acknowledgement would be appreciated.
%
% Extract all contours from CMR42 XML workspace document
% Uses java package XOM to parz file
% File parsing is fast, slow part is interpreting each element
% Structure of cvi42 XML file
% workspace - nodes / Workspace::MapHash
%   StudyMapStates - Workspace::MapList
%     StudyUid - QString
%     ImageStates - nodes / Workspace::MapHash
%       List all dicom images in study keyed by dicom SOP UID
%           a number of display parameters (ignored)
%           Contours - Workspace::MapHash - List of contours
%            List item n
%               name
%               ImageSize
%               Points               - actual contours
%               PixelSize
%               PointsBeforeRotation
%               SubpixelResolution   - scaling from points to dicom
%               Type - how was it drawn
%               IsManuallyDrawn
%               RotationAngle
%               Label                - user defined name
%     

% Create XOM objects
    [pathstr,~,~] =  fileparts(mfilename('fullpath'));
    javaaddpath([pathstr filesep 'XOM' filesep 'xom-1.2.10.jar']);
    import nu.xom.*
    import java.io.*
    
    
    %format_oct2019 = true;  % change to false if pre oct 2019
    try
    jfn = java.io.FileInputStream(fname);
    catch
    error('!!! ATTENTION !!! AN ERROR OCCURED: Please check the contour name is correct.');
    end

    % should stream it to remove \n and spaces
    parser = nu.xom.Builder(false);

    xmldoc = parser.build(jfn);

    % The first child should be workspace
    workspace = xmldoc.getChild(0);

    elist = workspace.getChildElements();

    % find StudyMapStates element
    for i = 0:elist.size()-1
        if strcmp(elist.get(i).getAttribute(0).getValue , 'StudyMapStates')
            StudyMapStates = elist.get(i);
            break
        end
    end

    if~exist('StudyMapStates','var')
        fprintf('\nCould not find StudyMapStates\n');
    end
    
    % get StudyUid and ImageStates node
    elist = StudyMapStates.getChild(1).getChildElements();
    
    for i = 0:elist.size()-1
        if strcmp(elist.get(i).getAttribute(0).getValue , 'ImageStates')
            ImageStates = elist.get(i);
            
        elseif strcmp(elist.get(i).getAttribute(0).getValue , 'StudyUid')
            nodeStudyUid = elist.get(i);
        end
    end
    
    if~exist('ImageStates','var') || ~exist('nodeStudyUid','var')
        fprintf('\nCould not find ImageStates\n');
    end
    study_uid = nodeStudyUid.getValue().toCharArray()';
    
    %For each UIDs listed extract its contour information from 'contours'
    % Can limit to getChildElements('item','http://www.circlecvi.com/cvi42/Workspace/Hash/')
    % but this is everything, so dont worry
    elist = ImageStates.getChildElements();   % List of all dicom images
    nImages = elist.size();

    contour_count = 0;
    for i = 0:nImages-1   % For each dicom UID
        elist_sub = elist.get(i).getChildElements();
%         if(format_oct2019)
%             elist_cont = elist_sub.get(2).getChildElements();  % 30 oct was ind 1
%         else
%             elist_cont = elist_sub.get(1).getChildElements();  % 30 oct was ind 1
%         end
        for i_sub = 0:elist_sub.size-1
            if(strcmp(elist_sub.get(i_sub).getAttribute(0).getValue(),'Contours'))
                elist_cont = elist_sub.get(i_sub).getChildElements();
                if(elist_cont.size() >0)
                    contour_count = contour_count+1;
                    contours(contour_count).uid  = elist.get(i).getAttribute(0).getValue().toCharArray()';
                    contours(contour_count).count = elist_cont.size();

                    % For each contour in this image
                    for i_cont = 0:elist_cont.size()-1
                        % Get ImageSize, Points, PixelSize, PointsBeforeRotation,
                        % SubpixelResolution, Type, IsManuallyDrawn, RotationAngle, Label
                        % loop through all points
                        el_contour = elist_cont.get(i_cont).getChildElements();
                        contours(contour_count).str(i_cont+1).name = elist_cont.get(i_cont).getAttribute(0).getValue().toCharArray()';
                        % For each 
                        for i_data = 0:el_contour.size()-1
                            [nodeVal, name] = getNode(el_contour.get(i_data));
                            contours(contour_count).str(i_cont+1).(name) = nodeVal;
                        end
            %             contours(contour_count).str(i_cont+1) = str;

                        % record contour name
                    end
                end
            end
        end
    end


end

function [nodeVal, name] = getNode(node)
    name = node.getAttribute(0).getValue().toCharArray()';
    
    switch(node.getAttribute(1).getValue().toCharArray()')
        case 'QSize'
            nodeVal = GetQSize(node);
        case 'QSizeF'
            nodeVal = GetQSize(node);
        case 'QPolygon'
            nodeVal = GetQPolygon(node);
        case 'qint16'
            nodeVal = GetDouble(node);
        case 'double'
            nodeVal = GetDouble(node);
        case 'bool'
            nodeVal = GetBool(node);
        case 'QString'
            nodeVal = GetQString(node);
        case 'Workspace::Cvi42::ContourType'
            nodeVal = GetQString(node);
        otherwise
    end

end

function sz = GetQSize(node)
    sz = [0 0];
    num_list = node.getValue().split('\n');
    sz(1) = str2double(num_list(2));
    sz(2) = str2double(num_list(3));
%     elist = node.getChildElements();
%     sz(1) = str2double( elist.get(0).getValue());
%     sz(2) = str2double( elist.get(1).getValue());
end

function poly = GetQPolygon(node)
% Use a split string as it is twice as fast as node traversal
    num_list = node.getValue().split('\n');
    poly = NaN(num_list.size());
    % order is space space number number space space
    for i = 3:4:num_list.size()-1
        poly(i-2) = str2double(num_list(i));
        poly(i-1) = str2double(num_list(i+1));
    end
    if(length(poly)>5)
        %minimum 6 elements for interpretation this way
        poly(3:4:end)=[];
        poly(3:3:end)=[];
        poly(end) = [];
        poly(end) = [];
        poly = reshape(poly,2,length(poly)/2);
    else
        poly = [];  % was an empty one
    end
%     chel = node.getChildElements();
%     poly = zeros(2,chel.size());
%     for i_pt = 1:chel.size()
%         poly(1,i_pt) = str2double(chel.get(i_pt-1).getChildElements().get(0).getValue());
%         poly(2,i_pt) = str2double(chel.get(i_pt-1).getChildElements().get(1).getValue());
%     end
end

function num = GetDouble(node)
    num = str2double(node.getValue());
end

function str = GetQString(node)
    str = node.getValue().toCharArray()';
end

function isTrue = GetBool(node)
    isTrue = false;    
    if(strcmp(node.getValue(),'true'))
        isTrue = true;
    end
end
