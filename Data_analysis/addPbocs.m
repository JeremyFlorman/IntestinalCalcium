fp = 'Z:\OAS\foodEncounter\noFood-0min';
%% pBoc
sheets = dir([fp '\**\*(dmp*.csv']);
if ~isempty(sheets)
    for i = 1:length(sheets)
        data = readtable(fullfile(sheets(i).folder,sheets(i).name));
        events = data{:,'Name'};
        frames = cellfun(@str2double, regexp(events,"\d*(?=-)","match", 'once'));

        wormDataFile = dir([sheets(i).folder '\*wormdata.mat']);

        load(fullfile(wormDataFile(1).folder, wormDataFile(1).name))
        wormdata.pBoc = frames;

        save(fullfile(wormDataFile(1).folder, wormDataFile(1).name), 'wormdata');

    end
end

%% exp
sheets = dir([fp '\**\*(exp*.csv']);
if ~isempty(sheets)
    for i = 1:length(sheets)
        data = readtable(fullfile(sheets(i).folder,sheets(i).name));
        events = data{:,'Name'};
        frames = cellfun(@str2double, regexp(events,"\d*(?=-)","match", 'once'));

        wormDataFile = dir([sheets(i).folder '\*wormdata.mat']);

        load(fullfile(wormDataFile(1).folder, wormDataFile(1).name))
        wormdata.exp = frames;

        save(fullfile(wormDataFile(1).folder, wormDataFile(1).name), 'wormdata');

    end
end

sheets = dir([fp '\**\*pumpingRate*.csv']);
if ~isempty(sheets)
    for i = 1:length(sheets)
        data = readtable(fullfile(sheets(i).folder,sheets(i).name));
        events = data{:,'Name'};
        pumpFrames = cellfun(@str2double, regexp(events,"\d*(?=-)","match", 'once'));

        wormDataFile = dir([sheets(i).folder '\*wormdata.mat']);

        load(fullfile(wormDataFile(1).folder, wormDataFile(1).name))


        %% Instantaneous frequency trace
        fs = 30;
        numFrames = length(wormdata.bulkSignal);
        pumpingRate = zeros(numFrames,1);

        for k = fs:numFrames
            lo = k-fs;
            hi = k;
            pumpsInWindow = pumpFrames(pumpFrames>lo & pumpFrames < hi);
            if ~isempty(pumpsInWindow)
                pumpingRate(k) = numel(pumpsInWindow);
            end
        end
        


        
        smoothedPumping = smoothdata(pumpingRate, "movmean", 60);


        if i == 7
            pumpingRate(6000:end) = nan;
        end

        wormdata.pumpTimes = frames;
        wormdata.pumpingRate = pumpingRate;

        save(fullfile(wormDataFile(1).folder, wormDataFile(1).name), 'wormdata');
        figure()
        plot(smoothedPumping)
    end
end