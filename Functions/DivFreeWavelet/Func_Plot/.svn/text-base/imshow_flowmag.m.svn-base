function h=imshow_flowmag(imMag,vx,vy,vz,vMax,orient,z)
%
% h = imshow_flowmag( imMag, vx, vy, vz, vMax, orient, z)
%
% Generic plotting function for flow magnitude
% The directions are oriented like matrix, ie (i,j,k) = (x,y,z)
% but should not matter for flow magnitude
% Modified from Joseph Cheng's imshow3s
%
% Inputs:
%     imMag   -   Image magnitude/image Mask
%     vx      -   Points "down" in matrix
%     vy      -   Points "right" in matrix
%     vz      -   Points "out" in matrix
%     vMax    -   Maximum velocity
%     orient  -   Figure orientation, defaults to [1,2,3] (optional)
%     z       -   Slice position (optional)
%
% Outputs:
%     h       -   Figure handle
%
% Example:
%     imshow_flowmag( imMag, vx, vy, vz, pi, [-3,-2,1] );
%     flips 2nd and 3rd dimension and permutes input data to order [3,2,1]
%
% (c) Frank Ong 2013

ha = fighandle();

ha.figure = gcf;
clf; % clear current figure.

if (~isreal(sum(imMag(:)))  )
    error('Input must be a real 3D volume');
end

vMag = sqrt(vx.^2+vy.^2+vz.^2);

if (nargin == 6 && ~isempty(orient) && sum((orient~=[1,2,3])))
    perm = abs(orient);
    imMag = permute(imMag,perm);
    vMag = permute(vMag,perm);
    
    
    for i = 1:3
        if (orient(i)<0)
            imMag = flipdim(imMag,i);
            vMag = flipdim(vMag,i);
        end
    end
end

if (nargin == 6 && size(imMag,3)==1)
    z = 1;
    ha.singleSlice = 1;
end

if (nargin == 7 && ~isempty(z))
    imMag = imMag(:,:,z);
    vMag = vMag(:,:,z);
    
    z = 1;
    ha.singleSlice = 1;
else
    z = round(size(imMag,3)/2);
    ha.singleSlice = 0;
    if size(imMag,3)==1
        ha.singleSlice = 1;
    end
end

%% Initialize constrast maps
ims = imMag(:,:,z);
map = [min(imMag(:)) max(imMag(:))]*1.0;
mapdiff = diff(map);
if mapdiff==0
    map = [map(1)-0.5,map(1)+0.5];
end
ha.map = map;
ha.imThresh = map(1)+mapdiff*0.2;

imMask = ims > (ha.imThresh);

%% Initialize Velocities
vmap = [0, vMax];
vmapdiff = vMax;
ha.vmap = vmap;

vs = vMag(:,:,z);

%% Setup figure and plot
figure(ha.figure); hold off;axis off;
axes('position', [0.11 0.11 0.7 0.8]);

cmap = jet(512);

implot1 = ims;
implot1(implot1>ha.map(2)) = ha.map(2);
implot1(implot1<ha.map(1)) = ha.map(1);
implot1 = (implot1-ha.map(1))/(ha.map(2)-ha.map(1));
implot2 = implot1;
implot3 = implot1;

implot = zeros(size(implot1,1),size(implot1,2),3);

cInd = vs(imMask);
cInd = (cInd - vmap(1))/(vmap(2)-vmap(1));
cInd = round(cInd*512);
cInd(cInd<1) = 1;
cInd(cInd>512) = 512;
implot1(imMask) = cmap(cInd,1);
implot2(imMask) = cmap(cInd,2);
implot3(imMask) = cmap(cInd,3);
implot(:,:,1) = implot1;
implot(:,:,2) = implot2;
implot(:,:,3) = implot3;

ha.fh = imshow(implot);

screen_size = get(0, 'ScreenSize');
figure_position = get(ha.figure, 'Position');
set(ha.figure, 'Position', [figure_position(1) figure_position(2) screen_size(3)/2 screen_size(4)/2])
axis off;


ha.a = gca;
ha.z = z;
ha.map = map;
ha.busy = 0;

if ~(ha.singleSlice)
    ht = uicontrol(ha.figure,'Style','text',...
        'Units', 'Normalized',...
        'Position',[0.07 0.05 0.2 0.05],...
        'String',sprintf('Slice: %d',z),'FontSize',14);
    
    uicontrol(ha.figure,'Style', 'slider',...
        'SliderStep',[1/size(imMag,3) 1/size(imMag,3)],...
        'Min',1,'Max',size(imMag,3),'Value',z,...
        'Units', 'Normalized',...
        'Position', [0.07 0.00 0.2 0.05],...
        'Callback', {@int_slider_slice_callback,ha,ht,imMag,vMag});
    
    
    
    imMax = max(imMag(:));
    imMin = min(imMag(:));
    ht = uicontrol(ha.figure,'Style','text',...
        'Units', 'Normalized',...
        'Position',[0.29 0.05 0.2 0.05],...
        'String',sprintf('IM Thresh: %.1f',ha.imThresh),'FontSize',14);
    
    uicontrol(ha.figure,'Style', 'slider',...
        'SliderStep',[1/50 1/50],...
        'Min',imMin,'Max',imMax,'Value',ha.imThresh,...
        'Units', 'Normalized',...
        'Position', [0.29 0.00 0.2 0.05],...
        'Callback', {@int_slider_imthresh_callback,ha,ht,imMag,vMag});
    
    
    hmin = uicontrol(ha.figure,'Style','text',...
        'Units', 'Normalized',...
        'Position',[0.51 0.05 0.2 0.05],...
        'String',sprintf('Min Vel: %.2f',vmap(1)),'FontSize',14);
    
    uicontrol(ha.figure,'Style', 'slider',...
        'Min',0,'Max',vmap(2)+vmapdiff,'Value',vmap(1),...
        'Units', 'Normalized',...
        'Position', [0.51 0.00 0.2 0.05],...
        'Callback', {@int_slider_min_callback,ha,hmin,imMag,vMag});
    
    
    hmax = uicontrol(ha.figure,'Style','text',...
        'Units', 'Normalized',...
        'Position',[0.73 0.05 0.2 0.05],...
        'String',sprintf('Max Vel: %.2f',vmap(2)),'FontSize',14);
    
    uicontrol(ha.figure,'Style', 'slider',...
        'Min',0,'Max',vmap(2)+vmapdiff,'Value',vmap(2),...
        'Units', 'Normalized',...
        'Position', [0.73 0.00 0.2 0.05],...
        'Callback', {@int_slider_max_callback,ha,hmax,imMag,vMag});
    
end

%% Setup output
if (nargout==1)
    h = ha.figure;
end

function int_slider_slice_callback(hslider,event,ha,ht,imMag,vMag)
z = round(get(hslider,'Value'));
set(hslider,'Value',z);
set(ht,'String',sprintf('Slice: %d',z));
ha.z = z;
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vMag);
    ha.busy = 0;
end


function int_slider_imthresh_callback(hslider,event,ha,ht,imMag,vMag)
ha.imThresh = get(hslider,'Value');
set(hslider,'Value',ha.imThresh);
set(ht,'String',sprintf('IM Thresh: %.1f',ha.imThresh));
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vMag);
    ha.busy = 0;
end


function int_slider_max_callback(hslider,event,ha,ht,imMag,vMag)

m = get(hslider,'Value');
set(hslider,'Value',m);
set(ht,'String',sprintf('Max Vel: %.2f',m));
ha.vmap = [min(ha.vmap(1),m-0.1) m];
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vMag);
    ha.busy = 0;
end


function int_slider_min_callback(hslider,event,ha,ht,imMag,vMag)

m = get(hslider,'Value');
set(hslider,'Value',m);
set(ht,'String',sprintf('Min Vel: %.2f',m));
ha.vmap = [m  max(ha.vmap(2),m+0.1)];
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vMag);
    ha.busy = 0;
end



function imshow3_plotflow(ha,imMag,vMag)
z = ha.z;
vmap = ha.vmap;

ims = imMag(:,:,z);
imMask = ims>(ha.imThresh);
vs = vMag(:,:,z);

cmap = jet(512);
implot1 = imMag(:,:,z);
implot1(implot1>ha.map(2)) = ha.map(2);
implot1(implot1<ha.map(1)) = ha.map(1);
implot1 = (implot1-ha.map(1))/(ha.map(2)-ha.map(1));
implot2 = implot1;
implot3 = implot1;

implot = zeros(size(implot1,1),size(implot1,2),3);

cInd = vs(imMask);
cInd = (cInd - vmap(1))/(vmap(2)-vmap(1));
cInd = round(cInd*512);
cInd(cInd<1) = 1;
cInd(cInd>512) = 512;
implot1(imMask) = cmap(cInd,1);
implot2(imMask) = cmap(cInd,2);
implot3(imMask) = cmap(cInd,3);
implot(:,:,1) = implot1;
implot(:,:,2) = implot2;
implot(:,:,3) = implot3;

set(ha.fh,'CData',implot);
drawnow;
