
 filename = 'C:\Users\Jeremy\Dropbox\NS&B2025\ventral bias alberto.mp4';

[fp, name, ext] = fileparts(filename);
info = imfinfo(filename);

fr = 15;    % frame rate of original video 
playrate = 2; % multiplier for playback framerate

startidx = 1 %460;
endidx =450 % 3638 %startidx+(30*fr);




outpath = fullfile(fp, [name '_' num2str(floor(startidx/15)) '-'...
    num2str(floor(endidx/fr)) '.mp4']);

% info = imfinfo(filename)
vW = VideoWriter(outpath, 'MPEG-4');
vW.Quality = 100;
vW.FrameRate = fr*playrate;

open(vW) 
figure()


for i = startidx:endidx
     img = uint8(imread(filename,i));
      imshow(img);
     writeVideo(vW,img)
end

 close(vW)
disp('DONE!!!')