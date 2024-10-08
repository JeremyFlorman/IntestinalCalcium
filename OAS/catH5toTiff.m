%% this function converts H5 videos from OAS to side-by-side tiff files, 
% similar to optosplit recordings

fld = 'C:\Users\Jeremy\Desktop\foodEncounter'


imgDir = dir([fld '\**\*behavior\*.h5']);
imgDir = unique({imgDir.folder});

for j = 1:length(imgDir)

    bfH5 = dir([imgDir{j} '\*.h5']);

    for i = 1:length(bfH5)
        bfFile = fullfile(bfH5(i).folder, bfH5(i).name);
        gcFile = strrep(bfFile, 'behavior','gcamp');
        tempBF = h5read(bfFile, '/data');
        tempGC = h5read(gcFile, '/data');
        if i == 1
            bfimg = tempBF;
            gcimg = tempGC;
        elseif i>1
            bfimg = cat(3,bfimg, tempBF);
            gcimg = cat(3,gcimg, tempGC);
        end
    end

    minlen = min(size(bfimg,3),size(gcimg,3));

    mergedImage = cat(2,gcimg(:,:,1:minlen), bfimg(:,:,1:minlen));



    outputFileName = strrep(imgDir{j}, 'flircamera_behavior', 'mergedImage.tif');
    for k = 1:length(mergedImage)
    imwrite(mergedImage(:, :, k), outputFileName, 'WriteMode', 'append','Compression','none');
    end

end
