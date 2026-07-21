% fp = 'C:\Users\Jeremy\Dropbox\ins-3 intestinal calcium fig\movies\'
% 
% fn = {[fp '211105_zfis178_wildtype-30mM-TA_3.tif'], [fp '211108_zfis178_wildtype-control_2.tif'], ...
%     [fp '220524_zfis178_tyra-3-30mM-TA_3.tif'],[fp '211110_zfis178_tyra-3-30mM-TA_2.tif'], ...
%     [fp '220525_zfis178_itr-1-30mM-TA_4.tif'], [fp '240507_zfis178_itr-1-30mM-TA_4.tif'],[fp '240913_zfis178_wildtype-30mM-TA_3-noFood.tif']};
% % 
% for f = length(fn)
%     filename = fn{f};
    filename = "C:\Users\AlkemaLab\Desktop\HisCat\4Well_10mMHA\N2_acr-2_unc-47_tag-168_10mM-HA.tif";
    [fp, name, ext] = fileparts(filename);
    info = imfinfo(filename);

    fr = 1;    % frame rate of original video
    playrate = 50; % multiplier for playback framerate

    startidx = 1; %460;
    endidx = length(info);



    outpath = strrep(filename, '.tif',['_' num2str(playrate) 'x_Speed.mp4']);

    % info = imfinfo(filename)
    vW = VideoWriter(outpath, 'MPEG-4');
    vW.Quality = 100;
    vW.FrameRate = fr*playrate;

    open(vW)
    


    for i = startidx:endidx
        img = uint8(imread(filename,i));
%         imshow(img);
        writeVideo(vW,img)
    end

    close(vW)
    disp('DONE!!!')
% end