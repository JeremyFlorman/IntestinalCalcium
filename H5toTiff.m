d = dir('Y:\OAS\5-HT\wildtype+5HT\');
dirflag= [d(:).isdir];
d = d(dirflag);

for j = 3:length(d)
    folder = fullfile(d(j).folder,d(j).name)
% folder = 'Y:\OAS\5-HT\wildtype+5HT\231006_zfis178_wildtype+5HT_1';

settings = returnPlotSettings;
settings.OAS = 1;
settings.traceylimit = [0 20];
settings.peakthreshold = 4;
settings.axylimit = [0 45];
settings.trimExperimentLength =1;

wd = dir([folder '\*wormdata.mat']);
h5file = dir([folder '\*behavior\*.h5']);

[mtdata, ~] = processWormdata(fullfile(wd(1).folder, wd(1).name), settings);

%%
timePre = settings.framerate*30;

locs = mtdata.peakLoc(2:end);
ints = mtdata.peakIntervals;
img =  [];
tic

for i = 1:length(h5file)
    tempfile = fullfile(h5file(i).folder, h5file(i).name);
    tempimg = h5read(tempfile, '/data');
    if i == 1
        img = tempimg;
    elseif i>1
        img = cat(3,img, tempimg);
    end
end
toc
%%


for i = 1:length(locs)
    outputFileName = [strrep(wd.name, 'wormdata.mat', ['bfImage_' num2str(i) '_interval(' num2str(round(ints(i))) ').tiff'])];
    localPath = ['C:\tmp\' outputFileName]
    destinationPath = [folder '\' outputFileName];

    temploc = locs(i);
    for k = temploc-timePre:temploc
        imwrite(img(:, :, k), localPath, 'WriteMode', 'append','Compression','none');
    end

    [status, ~, ~] = copyfile(localPath, destinationPath);
    
    if status==1
        disp([outputFileName ' Uploaded Successfully, cleaning up...'])
        delete(localPath);
    end

end
toc
end
