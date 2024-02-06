function h=imshow_flow(imMag,vx,vy,vz,vMax,orient,z)
%
% h = imshow_flow ( imMag, vx, vy, vz, vMax, orient, z)
%
% Generic plotting function for flow fields
% The directions are oriented like matrix, ie (i,j,k) = (x,y,z)
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
%     imshow_flow( imMag, vx, vy, vz, pi, [-3,-2,1] );
%     flips 2nd and 3rd dimension and permutes input data to order [3,2,1]
%
% (c) Frank Ong 2013

ha = fighandle();

ha.figure = gcf;
clf; % clear current figure.

if (~isreal(sum(imMag(:)))  )
    error('Input must be a real 3D volume');
end

if (nargin >= 6 && ~isempty(orient) && sum((orient~=[1,2,3])))
    perm = abs(orient);
    imMag = permute(imMag,perm);
    vx = permute(vx,perm);
    vy = permute(vy,perm);
    vz = permute(vz,perm);
    
    vxOrig = vx;
    vyOrig = vy;
    vzOrig = vz;
    if (perm(1)==2)
        vx = vyOrig;
    end
    if (perm(1)==3)
        vx = vzOrig;
    end
    if (perm(2)==1)
        vy = vxOrig;
    end
    if (perm(2)==3)
        vy = vzOrig;
    end
    if (perm(3)==1)
        vz = vxOrig;
    end
    if (perm(3)==2)
        vz = vyOrig;
    end
    clear vxOrig vyOrig vzOrig
    
    for i = 1:3
        if (orient(i)<0)
            imMag = flipdim(imMag,i);
            vx = flipdim(vx,i)*(2*(i~=1)-1);
            vy = flipdim(vy,i)*(2*(i~=2)-1);
            vz = flipdim(vz,i)*(2*(i~=3)-1);
        end
    end
end

if (nargin == 7 && ~isempty(z))
    imMag = imMag(:,:,z);
    vx = vx(:,:,z);
    vy = vy(:,:,z);
    vz = vz(:,:,z);
    z = 1;
    ha.singleSlice = 1;
else
    z = round(size(imMag,3)/2);
    ha.singleSlice = 0;
    if size(imMag,3)==1
        ha.singleSlice = 1;
    end
end

ha.subsample = min(max(round(log2(max(size(imMag))))-4,1),4);

%% Initialize constrast maps
ims = imMag(:,:,z);
map = [min(imMag(:)) max(imMag(:))]*1.0;
mapdiff = diff(map);
if mapdiff==0
    map = [map(1)-0.5,map(1)+0.5];
end
ha.map = map;
ha.imThresh = map(1)+mapdiff*0.2;


%% Initialize Velocities

vmap = [0, vMax];
ha.vmap = vmap;

ha.scale_vec = 8*ha.subsample/vMax;
[xgrid,ygrid] = meshgrid(1:size(imMag,2),1:size(imMag,1));
ha.xgrid = xgrid(1:ha.subsample:end,1:ha.subsample:end);
ha.ygrid = ygrid(1:ha.subsample:end,1:ha.subsample:end);

imMask = ims(1:ha.subsample:end,1:ha.subsample:end) > (ha.imThresh);
vxs = vx(1:ha.subsample:end,1:ha.subsample:end,z);
vys = vy(1:ha.subsample:end,1:ha.subsample:end,z);
vzs = vz(1:ha.subsample:end,1:ha.subsample:end,z);
vs = sqrt(vxs.^2+vys.^2+vzs.^2);

%% Setup figure and plot
figure(ha.figure); hold off;axis off;
axes('position', [0.11 0.11 0.7 0.8]);

ha.fh = imshow(imMag(:,:,z),ha.map);

hold on;

cmap = jet(32);
crange = ((0:31)/31);
vs = vs/ha.vmap(2);
ha.fh_flow = cell(1,32);
for i = 1:32
    if i==1
        vMask = vs<=crange(1);
    elseif i==32
        vMask = vs>=crange(31);
    else
        vMask = (vs<=crange(i))&(vs>=crange(i-1));
    end
    mask = vMask&imMask;
    ha.fh_flow{i} = quiver(ha.xgrid(mask),ha.ygrid(mask),vys(mask)*ha.scale_vec,vxs(mask)*ha.scale_vec,'Color',cmap(i,:),'AutoScale','off');
end
hold off;

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
        'Callback', {@int_slider_slice_callback,ha,ht,imMag,vx,vy,vz});
    
    
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
        'Callback', {@int_slider_imthresh_callback,ha,ht,imMag,vx,vy,vz});
    
    
    hscale = uicontrol(ha.figure,'Style','text',...
        'Units', 'Normalized',...
        'Position',[0.51 0.05 0.2 0.05],...
        'String',sprintf('Scale Vec: %.2f',ha.scale_vec),'FontSize',14);
    
    uicontrol(ha.figure,'Style', 'slider',...
        'Min',0,'Max',ha.scale_vec*4,'Value',ha.scale_vec,...
        'Units', 'Normalized',...
        'Position', [0.51 0.00 0.2 0.05],...
        'Callback', {@int_slider_scale_callback,ha,hscale,imMag,vx,vy,vz});
    
    
    hsub = uicontrol(ha.figure,'Style','text',...
        'Units', 'Normalized',...
        'Position',[0.73 0.05 0.2 0.05],...
        'String',sprintf('Subsample Vec: %d',ha.subsample),'FontSize',14);
    
    uicontrol(ha.figure,'Style', 'slider',...
        'SliderStep',[1/9 1/9],...
        'Min',1,'Max',10,'Value',ha.subsample,...
        'Units', 'Normalized',...
        'Position', [0.73 0.00 0.2 0.05],...
        'Callback', {@int_slider_sub_callback,ha,hsub,imMag,vx,vy,vz});
end


%% Setup output
if (nargout==1)
    h = ha.figure;
end

function int_slider_slice_callback(hslider,event,ha,ht,imMag,vx,vy,vz)
z = round(get(hslider,'Value'));
set(ht,'String',sprintf('Slice: %d',z));
ha.z = z;
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vx,vy,vz);
    ha.busy = 0;
end



function int_slider_imthresh_callback(hslider,event,ha,ht,imMag,vx,vy,vz)
ha.imThresh = get(hslider,'Value');
set(ht,'String',sprintf('IM Thresh: %.1f',ha.imThresh));
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vx,vy,vz);
    ha.busy = 0;
end

function int_slider_scale_callback(hslider,event,ha,ht,imMag,vx,vy,vz)

ha.scale_vec = get(hslider,'Value');
set(ht,'String',sprintf('Scale Vec: %.2f',ha.scale_vec));
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vx,vy,vz);
    ha.busy = 0;
end


function int_slider_sub_callback(hslider,event,ha,ht,imMag,vx,vy,vz)

ha.subsample = round(get(hslider,'Value'));
[xgrid,ygrid] = meshgrid(1:size(imMag,2),1:size(imMag,1));
ha.xgrid = xgrid(1:ha.subsample:end,1:ha.subsample:end);
ha.ygrid = ygrid(1:ha.subsample:end,1:ha.subsample:end);
set(ht,'String',sprintf('Subsample Vec: %d',ha.subsample));
if ~ha.busy
    ha.busy = 1;
    imshow3_plotflow(ha,imMag,vx,vy,vz);
    ha.busy = 0;
end



function imshow3_plotflow(ha,imMag,vx,vy,vz)
z = ha.z;
ims = imMag(:,:,z);
imMask = ims(1:ha.subsample:end,1:ha.subsample:end)>(ha.imThresh);
vxs = vx(1:ha.subsample:end,1:ha.subsample:end,z);
vys = vy(1:ha.subsample:end,1:ha.subsample:end,z);
vzs = vz(1:ha.subsample:end,1:ha.subsample:end,z);
vs = sqrt(vxs.^2+vys.^2+vzs.^2);
vs = vs/ha.vmap(2);

%figure(ha.figure);
set(ha.fh,'CData',ims);
set(ha.a, 'CLim',ha.map);
crange = ((0:31)/31);
for i = 1:32
    if i==1
        vMask = vs<crange(1);
    elseif i==32
        vMask = vs>=crange(31);
    else
        vMask = (vs<crange(i))&(vs>=crange(i-1));
    end
    mask = vMask&imMask;
    set(ha.fh_flow{i},'XData',1e-6);
    set(ha.fh_flow{i},'YData',1e-6);
    set(ha.fh_flow{i},'UData',1e-6);
    set(ha.fh_flow{i},'VData',1e-6);
    set(ha.fh_flow{i},'XData',ha.xgrid(mask));
    set(ha.fh_flow{i},'YData',ha.ygrid(mask));
    set(ha.fh_flow{i},'UData',vys(mask)*ha.scale_vec);
    set(ha.fh_flow{i},'VData',vxs(mask)*ha.scale_vec);
end
drawnow;
