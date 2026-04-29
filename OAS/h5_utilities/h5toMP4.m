h5path = 'C:\Users\AlkemaLab\Desktop\HisCat\zfex1028_Tag-168HisCat\Movies\260427_zfex1028_tag-168HisCat_+30mM_HA_1\2026_04_27_15_24_28_flircamera_behavior';
d = dir([h5path '\*.h5']);

[indices] = indexCameraTimestamps(h5path, 0);

fr = 15;    % frame rate of original video
playrate = 10; % multiplier for playback framerate

startidx = 1; %460;
endidx = length(indices.timestamps);


outpath = strrep(indices.h5_path{1}, '.h5',['_' num2str(playrate) 'x_Speed.mp4']);

vW = VideoWriter(outpath, 'MPEG-4');
vW.Quality = 100;
vW.FrameRate = fr*playrate;

open(vW)



for i = startidx:endidx
    img = getSlice(i, indices);
    % imshow(img);
    % drawnow()
    writeVideo(vW,img)
end

close(vW)
disp('DONE!!!')
% end