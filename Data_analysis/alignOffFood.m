function [exitInt1Matrix, exitInt9Matrix] = alignOffFood(wormdata) %
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
    wormdata = evalin("caller",'wormdata');
end

foodExits = wormdata.offFood;
foodEntries = wormdata.onFood;
secPost = 90;
secPre = 30;
fps = 15;

minExitDuration = 45*fps; % only consider food exits that last 45s or longer
minSpikeInterval = 70*fps; % exclude exit events with a spike that occurs a minute later
maxGap = 300*fps; % allow spikes only up to 300 sec after exit

intStart = 20;
intSegmentSize = size(wormdata.autoAxialSignal, 2)*.1;

exitTickLabels = {};
exitTraceCount = 1;
exitTraceMatrix = [];
exitInt1Matrix = [];
exitInt9Matrix = [];

entryTickLabels = {};
entryTraceCount = 1;
entryTraceMatrix = [];
entryInt1Matrix = [];
entryInt9Matrix = [];




%% Get calcium traces when animals leave food and align to last Ca2+ spike. Exclude short forays. 
for i = 1:length(foodExits)
    exitFrame = foodExits(i);
    entryIdx = find(foodEntries>exitFrame, 1, "first");
    entryFrame = foodEntries(entryIdx);
    framesOffFood = entryFrame-exitFrame;

    % exclude short exits
    if framesOffFood >= minExitDuration

        validSpikes = wormdata.peakLoc(wormdata.peakLoc > exitFrame & wormdata.peakLoc < entryFrame);
        validSpikes = validSpikes(validSpikes - exitFrame < maxGap); % keep only spikes within maxGap after exit

        if isempty(validSpikes)
            [~, spikeIdx] = min(abs(wormdata.peakLoc - exitFrame)); %if there were no spikes off food, use the nearest one on food
        else
            lastSpike = validSpikes(end);
            [~, spikeIdx] = min(abs(wormdata.peakLoc - lastSpike));
        end

        lastSpike = wormdata.peakLoc(spikeIdx); % the frame number of the last valid spike off food


        if spikeIdx ~= length(wormdata.peakLoc)
            subsequentSpike = wormdata.peakLoc(spikeIdx+1);
            framesBetweenSpikes = subsequentSpike-lastSpike; % check to see when the next spike is
        else
            framesBetweenSpikes = inf; % if there are no more spikes give up
        end

        if framesBetweenSpikes > minSpikeInterval

            preIdx = lastSpike - secPre*fps;
            postIdx = lastSpike + secPost*fps;

            if preIdx <1
                preIdx = 1;
            end

            if postIdx > size(wormdata.autoAxialSignal,1)
                postIdx =size (wormdata.autoAxialSignal,1);
            end

            %% Get and concatenate valid traces
            backgroundTrace = wormdata.backgroundSignal(preIdx:postIdx, :);
            exitTrace = wormdata.autoAxialSignal(preIdx:postIdx, :)-backgroundTrace;
            buffer = nan(size(exitTrace,1),1);
            exitTraceMatrix = horzcat(exitTraceMatrix, buffer, exitTrace);
            %% Get int1 and int9 signal
            int1Signal = exitTrace(:, intStart:intStart+intSegmentSize);
            int9Signal = exitTrace(:, end-intSegmentSize:end);
            exitInt1Matrix = horzcat(exitInt1Matrix, mean(int1Signal,2, "omitmissing"));
            exitInt9Matrix = horzcat(exitInt9Matrix, mean(int9Signal,2, "omitmissing"));

            exitTickLabels(exitTraceCount) = {num2str(exitTraceCount)};
            exitTraceCount = exitTraceCount+1;
        end
    end
end

for i =1:length(foodEntries)
    entryFrame = foodEntries(i);

    if entryFrame>1 % ignore the entry event if it was the first frame, this is not always annotated 
        % Find how long the off-food bout was
        if i==1 && min(foodEntries)<min(foodExits) % if the starting condition wasn't annotated and the earliest event was an entry, assume the worm started off-food
            exitFrame = 1; 
        else
            exitIdx = find(foodExits<entryFrame, 1, 'last'); % otherwise just find the previous exit event
            exitFrame = foodExits(exitIdx);
        end

        framesOffFood = entryFrame-exitFrame;
        timeOffFood = framesOffFood/fps
        
        if framesOffFood >= 90*fps % ensure the time off food was of sufficient duration

            preIdx = entryFrame - secPre*fps;
            postIdx = entryFrame + secPost*fps;

            if preIdx <1
                preIdx = 1;
            end

            if postIdx > size(wormdata.autoAxialSignal,1)
                postIdx =size (wormdata.autoAxialSignal,1);
            end
            %% Get and concatenate valid traces
            backgroundTrace = wormdata.backgroundSignal(preIdx:postIdx, :);
            entryTrace = wormdata.autoAxialSignal(preIdx:postIdx, :)-backgroundTrace;
            buffer = nan(size(entryTrace,1),1);
            entryTraceMatrix = horzcat(entryTraceMatrix, buffer, entryTrace);
            %% Get int1 and int9 signal
            int1Signal = entryTrace(:, intStart:intStart+intSegmentSize);
            int9Signal = entryTrace(:, end-intSegmentSize:end);
            entryInt1Matrix = horzcat(entryInt1Matrix, mean(int1Signal,2, "omitmissing"));
            entryInt9Matrix = horzcat(entryInt9Matrix, mean(int9Signal,2, "omitmissing"));

            entryTickLabels(entryTraceCount) = {num2str(entryTraceCount)};
            entryTraceCount = entryTraceCount+1;

        end
    end
end


figure('Position', [347.4000 61 631.2000 683.2000], 'Color', [1 1 1]);
% t = tiledlayout(3,2);
%% Exit Traces

% exitAx = nexttile([3 1]);
traceX = linspace(-secPre, secPost, size(exitTraceMatrix,1));
traceY = 1:size(exitTraceMatrix,2)+(exitTraceCount-1);
imagesc(traceX, traceY, smoothdata(exitTraceMatrix, 'gaussian', 60)', [0 90])

yt = nan(exitTraceCount-1,1);
yt(1) = 100;
if exitTraceCount >1
    for i = 2:exitTraceCount-1
        yt(i) = yt(i-1)+201;
    end
end

yticks(yt);
yticklabels(exitTickLabels)
colormap("turbo");

ylabel('Leaving Event #')
xlabel('Time since last Ca^2^+ wave (s)')
xlim([-15 60])

xticks(-15:15:60)


%% Entry Traces

% entryAx = nexttile([3 1]);
% traceX = linspace(-secPre, secPost, size(entryTraceMatrix,1));
% traceY = 1:size(entryTraceMatrix,2)+(entryTraceCount-1);
% imagesc(traceX, traceY, smoothdata(entryTraceMatrix, 'gaussian', 60)', [0 90])
% 
% yt = nan(entryTraceCount-1,1);
% yt(1) = 100;
% if entryTraceCount >1
%     for i = 2:entryTraceCount-1
%         yt(i) = yt(i-1)+201;
%     end
% end
% 
% yticks(yt);
% yticklabels(entryTickLabels)
% colormap("turbo");
% 
% ylabel('Entry Event #')
% xlabel('Time since food entry (s)')
xlim([-15 60])
cb = colorbar;
cb.Label.String='GCaMP Signal (a.u.)';
cb.Label.Rotation = -90;


cb = colorbar;
cb.Label.String='GCaMP Signal (a.u.)';
cb.Label.Rotation = -90;