function [phaseData] = computePhaseChange(wormdata, fps, minOffDurationSeconds, h5path, dilationInMM, debugPlots)
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
    dilationInMM = -0.1;
    debugPlots = 1;
elseif nargin < 4
    dilationInMM = -0.1;
    debugPlots = 1;
end

%% Recompute Food Vector
% dilationInMM = -0.1;
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
        secondsRemaining(secondsRemaining<0) = 0;
        frameOfNextEvent  = loc(find(loc>frameOfThisEvent, 1)); % find the next ca2+ wave after the leaving event
        isNextEventOnFood = onFoodVector(frameOfNextEvent); % check the food vector to see if it was on or off food

        if isNextEventOnFood == 1
            eventBout = onBouts(find(onBouts(:,1)<frameOfNextEvent, 1, "last"),:); % get the bout where the next event happens
            entryFrame =eventBout(1);
            secondsAfterFoodEntry = (frameOfNextEvent-entryFrame)/fps; % find how many seconds after food entry the next event occured
            combinedCycleTime = secondsRemaining+secondsAfterFoodEntry; % cycle time excluding time off food

            phaseChange = secondsAfterFoodEntry - secondsRemaining; % this is signed, will tell you +/- cycle extension/reduction

            exitTrace = wormdata.autoAxialSignal(frameOfThisEvent-30*fps:frameOfThisEvent+90*fps,:); % Get signal when worm exits food, aligned to last Ca2+ wave
            [int1Exit, int9Exit] = extractCellSignal(exitTrace); % get int1 and int9 on exit

            entryTrace = wormdata.autoAxialSignal(entryFrame-30*fps:entryFrame+90*fps,:); % get signal when worm re-enters food
            [int1Entry, int9Entry] = extractCellSignal(entryTrace); % get int1 and int9 on reentry

            int1EntryWarped = interpolateKymograph(int1Entry, 30*fps, 30*fps+(frameOfNextEvent-entryFrame)); % interpolate int1/9 signal. traces begin 30sec before food entry so 30*fps is the first timepoint
            int9EntryWarped = interpolateKymograph(int9Entry, 30*fps, 30*fps+(frameOfNextEvent-entryFrame)); % (frameOfNextEvent-entryFrame) is delay to the first Ca2+ wave, so we add this to 30*fps to get next timepoint

            % Store data
            fn = strsplit(wormdata.filename, '\');
            phaseData(eventIndex).filename = fn{end}; % file name of recording
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
            
            % Axial and Int1/Int9 Signal 
            phaseData(eventIndex).exitTrace = {exitTrace};
            phaseData(eventIndex).entryTrace = {entryTrace};
            
            phaseData(eventIndex).int1Exit = {int1Exit};
            phaseData(eventIndex).int9Exit = {int9Exit};

            phaseData(eventIndex).int1Entry = {int1Entry};
            phaseData(eventIndex).int9Entry = {int9Entry};
            
            phaseData(eventIndex).int1EntryWarped = {int1EntryWarped};
            phaseData(eventIndex).int9EntryWarped = {int9EntryWarped};

            eventIndex = eventIndex+1;

        end
    end
end

%% Plot for debugging
if debugPlots == 1
    if ~isempty(h5path) && isfield(phaseData, 'exitFrame')
        figure('Name', phaseData(1).filename, 'Position', [285.8000 89.8000 917.6000 654.4000], 'Color', [1 1 1]);
        t = tiledlayout(3, length(phaseData)+2, 'TileSpacing','compact','Padding','tight');
        h5Files = dir([h5path '\*.h5']);
        %% Display images of leaving and return events
        for i = 1:length(phaseData)
            nexttile(i)
            exitImg = getImage(phaseData(i).exitFrame, h5Files);
            imshow(exitImg)
            line(wormdata.headLoc(phaseData(i).exitFrame,1), wormdata.headLoc(phaseData(i).exitFrame,2), 'Color', 'g', 'Marker', 'o')
            line(wormdata.tailLoc(phaseData(i).exitFrame,1), wormdata.tailLoc(phaseData(i).exitFrame,2), 'Color', 'r', 'Marker', 'o')
            if i == 1
                ylabel('Exit')
                yticks([])
            end
            nexttile(length(phaseData)+i+2)
            entryImg = getImage(phaseData(i).entryFrame, h5Files);
            imshow(entryImg)
            line(wormdata.headLoc(phaseData(i).entryFrame,1), wormdata.headLoc(phaseData(i).entryFrame,2), 'Color', 'g', 'Marker', 'o')
            line(wormdata.tailLoc(phaseData(i).entryFrame,1), wormdata.tailLoc(phaseData(i).entryFrame,2), 'Color', 'r', 'Marker', 'o')
            if i == 1
                ylabel('Return')
                yticks([])
            end

        end
        %% plot bulk signal and event locations
        nexttile(((length(phaseData)+2)*2)+1, [1 length(phaseData)+2])
        patchAlpha = 1;
        patchColor = [0.996 0.9400 0.7920]; %[0.93 0.69 0.13]
        pk = max(bulkSignal, [], "all");
        if isfield(phaseData, 'preFrame')
            preSpikes = time(vertcat(phaseData.preFrame));
            postSpikes = time(vertcat(phaseData.postFrame));
            plot(time,bulkSignal, 'k',preSpikes,pk*1.05, 'rv',postSpikes,pk*1.05, 'gv', 'MarkerSize',3)
            text()
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

    %% Plot Exit and entry traces
    exTraces = [];
    enTraces = [];
    if isfield(phaseData, 'exitTrace')
        for i =1:numel(phaseData)
            buffer = nan(1,size(phaseData(i).exitTrace{:}',2));
            exTraces = vertcat(exTraces,buffer, phaseData(i).exitTrace{:}');
            enTraces = vertcat(enTraces,buffer, phaseData(i).entryTrace{:}');
        end

        traceX = linspace(-30, 90, size(exTraces,1));
        traceY = 1:size(exTraces,1)+(eventIndex-1);

        yt = nan(eventIndex-1,1);
        yt(1) = 100;
        yLabels = '1';
        if eventIndex >1
            for i = 2:eventIndex-1
                yt(i) = yt(i-1)+201;
                yLabels = [yLabels; num2str(i)];
            end
        end


        % figure();
        % tiledlayout(2,2, "TileSpacing","compact", 'Padding','tight');
        nexttile([2 1])
        imagesc(traceX, traceY, smoothdata(exTraces,2,'gaussian', 60), [0 45]);
        yticks(yt);
        yticklabels(yLabels)
        xticks(-30:30:90)
        box off

        nexttile([2 1])
        imagesc(traceX, traceY, smoothdata(enTraces,2,'gaussian', 60), [0 45])
        yticks(yt);
        ax = gca; 
        ax.YAxis.Visible = 'off';
        xticks(-30:30:90)
        
        box off
        % colormap("turbo")
    end


end




end

function img = getImage(frame, h5Files)
fileIndex  = floor((frame-1)/3600) + 1;
frameIndex = mod(frame-1, 3600) + 1;
filepath = fullfile(h5Files(fileIndex).folder, h5Files(fileIndex).name);
img = h5read(filepath, '/data',[1 1 frameIndex],[512 512 1]);
end