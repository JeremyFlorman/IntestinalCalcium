function [] = makeVideoFromTrack(path, timepoint, window)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
% 
% 
% % Initialize data cursor object
% cursorobj = datacursormode(gcf);
% cursorobj.SnapToDataVertex = 'on'; % Snap to our plotted data, on by default
% 
% while ~waitforbuttonpress 
%     % waitforbuttonpress returns 0 with click, 1 with key press
%     % Does not trigger on ctrl, shift, alt, caps lock, num lock, or scroll lock
%     cursorobj.Enable = 'on'; % Turn on the data cursor, hold alt to select multiple points
% end
% cursorobj.Enable = 'off';
% 
% mypoints = getCursorInfo(cursorobj);
% timepoint = mypoints(3);
enableDefaultInteractivity(gca);

[xPt,yPt] = ginput(1);
% Get (x,y) coordinates for all points
h = gco();
hx = h.XData;
hy = h.YData;
hz = h.ZData;
% Find the nearest point to selection
d = sqrt((xPt-hx).^2 + (yPt-hy).^2);
[~,minIdx] = min(d);

timepoint = hz(minIdx);



txt = input(['Save video at time: ' num2str(timepoint/15/60) ' min? (y/n)...'],"s");
if strcmp(txt,'y')
    makeVideoFromTimepoints(path, timepoint, window)
end
end