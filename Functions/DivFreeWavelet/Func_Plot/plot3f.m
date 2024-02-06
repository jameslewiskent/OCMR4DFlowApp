function h=plot3f(im,orient,z)
%
% h = plot3f( imMag, orient, z)
%
% Plotting function that plots velocity profiles with x and z slice control
% Modified from Joseph Cheng's imshow3s
% 
% Inputs: 
%     im      -   3D data
%     orient  -   Figure orientation, defaults to [1,2,3] (optional) 
%     z       -   Z slice position (optional) 
%
% Outputs:
%     h       -   Figure handle
%       
% (c) Frank Ong 2013

ha = fighandle();

ha.figure = gcf;
clf; % clear current figure.

if (~isreal(sum(im(:)))  )
    im = abs(im);
end

if (nargin >= 2 && ~isempty(orient) && sum((orient~=[1,2,3])))
    perm = 1:length(size(im));
    perm(1:length(orient)) = abs(orient);
    im = permute(im,perm);
    
    for i = 1:length(orient)
        if (orient(i)<0)
            im = flipdim(im,i);
        end
    end
end

if (nargin == 3 && ~isempty(z))
    im = im(:,:,z);
    z = 1;
    ha.singleSlice = 1;
else
    z = round(size(im,3)/2);
    ha.singleSlice = 0;
    if (size(im,3)==1)
        ha.singleSlice=1;
    end
end

x = round(size(im,1)/2);

ps = squeeze(im(x,:,z));

%% Setup figure and plot
figure(ha.figure); hold off;axis off;
axes('position', [0.11 0.11 0.7 0.8]);

ha.fh = plot(ps);
axis([1 length(ps) min(im(:)) max(im(:))]);

screen_size = get(0, 'ScreenSize');
figure_position = get(ha.figure, 'Position');
set(ha.figure, 'Position', [figure_position(1) figure_position(2) screen_size(3)/2 screen_size(4)/2])

ha.a = gca;
ha.z = z;
ha.x = x;


if ~(ha.singleSlice)
    plot3(ha,im);
    hslice = uicontrol(ha.figure,'Style','text',...
        'Units', 'Normalized',...
        'Position',[0.6 0.05 0.2 0.05],...
        'String',sprintf('Z Slice: %d',z),'FontSize',14);
    
    uicontrol(ha.figure,'Style', 'slider',...
        'SliderStep',[1/size(im,3) 1/size(im,3)],...
        'Min',1,'Max',size(im,3),'Value',z,...
        'Units', 'Normalized',...
        'Position', [0.6 0.00 0.2 0.05],...
        'Callback', {@int_slider_slice_callback,ha,hslice,im});
end


hx = uicontrol(ha.figure,'Style','text',...
    'Units', 'Normalized',...
    'Position',[0.20 0.05 0.2 0.05],...
    'String',sprintf('X Slice: %d',x),'FontSize',14);

uicontrol(ha.figure,'Style', 'slider',...
        'Min',1,'Max',size(im,1),'Value',x,...
    'Units', 'Normalized',...
    'Position', [0.20 0.00 0.2 0.05],...
    'Callback', {@int_slider_x_callback,ha,hx,im});



%% Setup output
if (nargout==1)
    h = ha.figure;
end



function int_slider_x_callback(hslider,event,ha,ht,im)

ha.x = round(get(hslider,'Value'));
set(hslider,'Value',ha.x);
set(ht,'String',sprintf('X Slice: %d',ha.x));
plot3(ha,im);


function int_slider_slice_callback(hslider,event,ha,ht,im)

z = round(get(hslider,'Value'));
set(hslider,'Value',z);
set(ht,'String',sprintf('Z Slice: %d',z));
ha.z = z;
plot3(ha,im);



function plot3(ha,im)

z = ha.z;
x = ha.x;

ps = squeeze(im(x,:,z));

set(ha.fh,'YData',ps);
drawnow;

