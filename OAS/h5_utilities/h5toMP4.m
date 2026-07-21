h5path = 'C:\src\OpenAutoScope-v2\data\260616_wildtype-10mM-HA\2026_06_16_14_55_43_flircamera_behavior';
d = dir([h5path '\*.h5']);

[indices] = indexCameraTimestamps(h5path, 0);

fr = 15;    % frame rate of original video
playrate = 5; % multiplier for playback framerate

startidx = 1; %460;
endidx = length(indices.timestamps);


outpath = strrep(indices.h5_path{1}, '.h5',['_' num2str(playrate) 'x_Speed.mp4']);

vW = VideoWriter(outpath, 'MPEG-4');
vW.Quality = 100;
vW.FrameRate = fr*playrate;

open(vW)

addTimestamp =1; 



%%
t = linspace(0, endidx/15/60, endidx);
for i = startidx:endidx
    if addTimestamp == 0
        img = getSlice(i, indices);
        writeVideo(vW,img)
    else
        img = getSlice(i, indices);
        imshow(img);
        txtstr = [num2str(t(i),'%0.2f') ' min'];
        text(5,7, txtstr, "FontSize", 12, "Color",[0 0 0])
        frame = getframe;
        % drawnow()
        writeVideo(vW,frame)
    end
end
% close(fig)
close(vW)
disp('DONE!!!')
% end