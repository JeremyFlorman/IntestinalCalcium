function [phaseData] = computePhaseChange(wormdata, fps, minOffDurationSeconds, h5path)
%computePhaseChange calculates the change in phase of the defecation cycle
%when a worm leaves and returns to food
%   inputs:
%   wormdata - wormdata structure containing bulkSignal and food
%   boutData
%
%   framerate (optional): frames per second of recording. Can be left
%   blank, or can be a scalar value indicating framerate. Or if plot
%   settings are passed, the framerate will be extracted from the settings
%   structure.
%



% Get Framerate from settings
if nargin<1
    wormdata = evalin("caller",'wormdata');
    fps = 15;
elseif nargin == 1
    fps = 15;
elseif nargin == 2 && isstruct(fps) % if we are passing the settings structure, find framerate
    fps = fps.framerate;
    minOffDurationSeconds = 15;
    h5path = [];
elseif nargin == 3
    h5path = [];
end

%% Recompute Food Vector
dilationInMM = -0.25;
onFoodVector = computeFoodVector(wormdata.patchROIs, dilationInMM);
%% Recalculate Food Bouts
boutData = computeFoodBouts(onFoodVector, fps, minOffDurationSeconds);

onBouts = boutData.onFood;
bulkSignal = wormdata.bulkSignal;% - wormdata.backgroundSignal;
loc = wormdata.peakLoc;
time = linspace(0, length(bulkSignal)/fps/60, length(bulkSignal));
phaseData = struct();

%% Remove short exits
% minOffDuration = 15;
%
% for i = 1:size(offBouts,1)
%     thisOff = offBouts(i,:);
%     secondsOff = (thisOff(2)-thisOff(1))/fps;
%     if secondsOff < minOffDuration
%         previousOnIdx = find(onBouts(:,2) == thisOff(1)); % find the previous bout based on the shared end/start times
%         nextOnIdx = find(onBouts(:,1) == thisOff(2)); % find the next bout based on the shared start/end times (I know, could just do +1 but I'm scared of edge cases)
%
%         onBouts(previousOnIdx, 2) = onBouts(nextOnIdx,2); % replace the end sof the previous onBout with the end of the next bout
%
%         offBouts(i,:) = [nan nan]; % clear the intervening offBout
%         onBouts(nextOnIdx,:) = [nan nan];  % clear the next onBout as it has been merged with the previous
%     end
% end



onFoodIntervals = [];
intervalIndex = 1;
%% Get Average Interval On-Food
for i = 1:size(onBouts,1)
    thisOn = onBouts(i,:);
    validEvents = loc(loc>thisOn(1) & loc<thisOn(2)); % find defecation events that occur during the current bout
    if numel(validEvents)>1
        onFoodIntervals(intervalIndex,1) = mean(diff(validEvents))/fps; % average interval during this bout
        intervalIndex = intervalIndex+1;
    end

end


%% Get Relative spike times
eventIndex = 1;
for i = 1:size(onBouts,1)
    thisOn = onBouts(i,:);
    exitFrame = thisOn(2);
    validEvents = loc(loc>thisOn(1) & loc<thisOn(2)); % find defecation events that occur during the current bout
    if numel(validEvents)>0

        interval = mean(onFoodIntervals, 1); % average interval during on-food bouts
        frameOfThisEvent = validEvents(end); % Frame of the last event during this on food bout
        secondsRemaining = interval - (exitFrame-frameOfThisEvent)/fps; % time between ca2+ wave and leaving event
        frameOfNextEvent  = loc(find(loc>frameOfThisEvent, 1)); % find the next ca2+ wave after the leaving event
        isNextEventOnFood = onFoodVector(frameOfNextEvent); % check the food vector to see if it was on or off food

        if isNextEventOnFood == 1
            eventBout = onBouts(find(onBouts(:,1)<frameOfNextEvent, 1, "last"),:); % get the bout where the next event happens
            entryFrame =eventBout(1);
            secondsAfterFoodEntry = (frameOfNextEvent-entryFrame)/fps; % find how many seconds after food entry the next event occured
            combinedCycleTime = secondsRemaining+secondsAfterFoodEntry; % cycle time excluding time off food

            phaseChange = combinedCycleTime -interval % this is signed, will tell you +/- cycle extension/reduction

            % Store data
            phaseData(eventIndex).preFrame = frameOfThisEvent; % Frame when the on-food defecation occured
            phaseData(eventIndex).exitFrame = exitFrame; % Frame when worm leaves food
            phaseData(eventIndex).entryFrame = entryFrame; % frame when worm returns to food
            phaseData(eventIndex).postFrame = frameOfNextEvent; % Frame of the subsequent defecation event after returning to food
            phaseData(eventIndex).interval = interval; % Average interval during the current onFood bout
            phaseData(eventIndex).secondsEventToExit = (exitFrame-frameOfThisEvent)/fps; % How many seconds between the on-food defecation event and the leaving event
            phaseData(eventIndex).totalCycleTime = (frameOfNextEvent-frameOfThisEvent)/fps; % Total time in seconds between the two defecation events
            phaseData(eventIndex).timeOffFood = (entryFrame - exitFrame)/fps; % Total time in seconds spent off food, this may also include on short on-food bouts where no defecation occurs
            phaseData(eventIndex).cycleSecondsRemaining = secondsRemaining; % The number of seconds left in the cycle (based on the average interval) when the animal left food
            phaseData(eventIndex).secondsDelay = secondsAfterFoodEntry; % The delay in seconds from food entry to the next defecation event
            phaseData(eventIndex).phaseChange = phaseChange; % The number of seconds the cycle shifted, simply (# seconds left in the cycle + delay after food entry) - average interval;
            eventIndex = eventIndex+1;

        end
    end
end

phaseTable = struct2table(phaseData);


%% Plot for debugging

if ~isempty(h5path)
    figure;
    t = tiledlayout(length(phaseData)+1, 2, 'TileSpacing','none','Padding','tight');
    h5Files = dir([h5path '\*.h5']);
    for i = 1:length(phaseData)
        nexttile()
        exitImg = getImage(phaseData(i).exitFrame, h5Files);
        imshow(exitImg)
        line(wormdata.headLoc(phaseData(i).exitFrame,1), wormdata.headLoc(phaseData(i).exitFrame,2), 'Color', 'g', 'Marker', 'o')
        line(wormdata.tailLoc(phaseData(i).exitFrame,1), wormdata.tailLoc(phaseData(i).exitFrame,2), 'Color', 'r', 'Marker', 'o')

        nexttile()
        entryImg = getImage(phaseData(i).entryFrame, h5Files);
        imshow(entryImg)
        line(wormdata.headLoc(phaseData(i).entryFrame,1), wormdata.headLoc(phaseData(i).entryFrame,2), 'Color', 'g', 'Marker', 'o')
        line(wormdata.tailLoc(phaseData(i).entryFrame,1), wormdata.tailLoc(phaseData(i).entryFrame,2), 'Color', 'r', 'Marker', 'o')
    end
    %% plot bulk signal and event locations
    nexttile([1 2])
    patchAlpha = 1;
    patchColor = [0.996 0.9400 0.7920]; %[0.93 0.69 0.13]
    pk = max(bulkSignal, [], "all");
    if isfield(phaseData, 'preFrame')
        preSpikes = time(vertcat(phaseData.preFrame));
        postSpikes = time(vertcat(phaseData.postFrame));
        plot(time,bulkSignal, 'k',preSpikes,pk*1.05, 'rv',postSpikes,pk*1.05, 'gv', 'MarkerSize',3)
    else
        plot(time,bulkSignal, 'k')
    end
    xpatch = nan(4,size(onBouts,1));
    for j=1:size(onBouts,1)
        s = onBouts(j,1);
        e = onBouts(j,2);
        xpatch(1:4, j) = [s s e e];
    end
    ypatch = repmat([0; 1; 1; 0],1 ,size(xpatch, 2));

    if isfield(wormdata, 'boutData') && ~isempty(wormdata.boutData)
        xCoords = xpatch/fps/60;
        yCoords = (ypatch*max(bulkSignal, [], 'all')*1.1)+0.1;
        p = patch(xCoords, yCoords, patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
        uistack(p, 'bottom');
    end


else
    %% Just plot bulk signal and event locations
    patchAlpha = 1;
    patchColor = [0.996 0.9400 0.7920]; %[0.93 0.69 0.13]
    pk = max(bulkSignal, [], "all");
    if isfield(phaseData, 'preFrame')
        preSpikes = time(vertcat(phaseData.preFrame));
        postSpikes = time(vertcat(phaseData.postFrame));
        plot(time,bulkSignal, 'k',preSpikes,pk*1.05, 'rv',postSpikes,pk*1.05, 'gv', 'MarkerSize',3)
    else
        plot(time,bulkSignal, 'k')
    end
    xpatch = nan(4,size(onBouts,1));
    for j=1:size(onBouts,1)
        s = onBouts(j,1);
        e = onBouts(j,2);
        xpatch(1:4, j) = [s s e e];
    end
    ypatch = repmat([0; 1; 1; 0],1 ,size(xpatch, 2));

    if isfield(wormdata, 'boutData') && ~isempty(wormdata.boutData)
        xCoords = xpatch/fps/60;
        yCoords = (ypatch*max(bulkSignal, [], 'all')*1.1)+0.1;
        p = patch(xCoords, yCoords, patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
        uistack(p, 'bottom');
    end
end


end

function img = getImage(frame, h5Files)
fileIndex  = floor((frame-1)/3600) + 1;
frameIndex = mod(frame-1, 3600) + 1;
filepath = fullfile(h5Files(fileIndex).folder, h5Files(fileIndex).name);
img = h5read(filepath, '/data',[1 1 frameIndex],[512 512 1]);
end