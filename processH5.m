function [bf, gfp, time] = processH5(foldername)
%Interpolates OAS behavior and GCaMP files to synchronize time
%   Detailed explanation goes here
% foldername = 'C:\src\OpenAutoScope-v2\data\zfis178';
d = dir([foldername '\*.h5']);

for i = 1:length(d)
    bPath = fullfile(d(i).folder,d(i).name);
    gPath =strrep(bPath,'behavior', 'gcamp');

    bTimes = h5read(bPath, '/times');
    gTimes = h5read(gPath, '/times');

    bData = h5read(bPath, '/data');
    gData = h5read(gPath, '/data');
    bFrames = size(bData,3);
    gFrames = size(gData,3);


    idx = NaN(gFrames,1);


    for j=1:bFrames
        ii = find(gTimes >=bTimes(j), 1);
        if ~isempty(ii)
            idx(j) = ii;
        end
    end
    % tic
    idx = idx(~isnan(idx));
    tempbf = bData(:,:,1:length(idx));
    tempgfp = gData(:,:,idx);
    temptime = bTimes(1:length(idx));
    if i == 1
        starttime = temptime(1);
    end
    temptime = temptime - starttime;

%     for j = 1:length(gData)
%         imshowpair(bf(:,:,j) ,imadjust(gfp(:,:,j) , [0.0 0.01]))
%         text(20,20, ['Time: ' num2str(round(time(j),2)) ' sec'])
%         pause(0.0001)
%     
%     end

    if i == 1
        bf = tempbf;
        gfp = tempgfp;
        time = temptime;
    elseif i>1
        bf = cat(3,bf,tempbf);
        gfp = cat(3,gfp, tempgfp);
        time = cat(1,time, temptime);
    end
    plot(time)
end




end