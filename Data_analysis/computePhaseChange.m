function [phaseChanges] = computePhaseChange(wormdata, framerate)
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
elseif isscalar(framerate)
    fps = framerate;
elseif isstruct(framerate)
    fps = framerate.framerate;
end


onBouts = wormdata.boutData.onFood;
offBouts = wormdata.boutData.offFood;
bulkSignal = wormdata.bulkSignal - wormdata.backgroundSignal;
loc = wormdata.peakLoc;
time = linspace(0, length(bulkSignal)/fps/60, length(bulkSignal));

%% Remove short exits
minOffDuration = 15;

for i = 1:size(offBouts,1) 
    thisOff = offBouts(i,:);
    secondsOff = (thisOff(2)-thisOff(1))/fps;
    if secondsOff < minOffDuration
        previousOnIdx = find(onBouts(:,2) == thisOff(1)); % find the previous bout based on the shared end/start times
        nextOnIdx = find(onBouts(:,1) == thisOff(2)); % find the previous bout based on the shared start/end times (I know, could just do +1 but I'm scared of edge cases)

        onBouts(previousOnIdx, 2) = onBouts(nextOnIdx,2); % replace the end of the previous onBout with the end of the next bout
        
        offBouts(i,:) = [nan nan]; % clear the intervening offBout
        onBouts(nextOnIdx,:) = [nan nan];  % clear the next onBout as it has been merged with the previous
    end
end


%% Get Relative spike times
for i = 1:size(onBouts,1)
    thisOn = onBouts(i,:);
    
    validEvents = loc(loc>thisOn(1) & loc<thisOn(2)); % find defecation events that occur during the current bout
    if numel(validEvents)>1
       
        interval = mean(diff(validEvents))/fps; % average interval during this bout
        frameOfThisEvent = validEvents(end); % Frame of the last event during this on food bout
        secondsRemaining = interval - (thisOn(2)-frameOfThisEvent)/fps; % time between ca2+ wave and leaving event
        frameOfNextEvent  = loc(find(loc>frameOfThisEvent, 1)); % find the next ca2+ wave after the leaving event
        isEventOnFood = wormdata.onFoodVector(frameOfNextEvent); % check the food vector to see if it was on or off food

        if isEventOnFood == 1
            eventBout = onBouts(find(onBouts(:,1)<frameOfNextEvent, 1, "last"),:); % get the bout where the next event happens
           
            secondsAfterFoodEntry = (frameOfNextEvent-eventBout(1))/fps; % find how many seconds after food entry the next event occured
           
            combinedCycleTime = secondsRemaining+secondsAfterFoodEntry; % cycle time excluding time off food
          
            % phaseChange = mod(combinedCycleTime,interval); %this is unsigned
           
            phaseChange = combinedCycleTime -interval % this is signed, will tell you +/- cycle extension
        end

        rawInterval = frameOfNextEvent-validEvents(end)/fps;
        rawPhaseChange = mod(rawInterval, interval)



        %% Plot for debugging
        patchAlpha = 1;
        patchColor = [0.996 0.9400 0.7920]; %[0.93 0.69 0.13]
        figure;
        pk = max(bulkSignal, [], "all");
        plot(time,bulkSignal, 'k',time([validEvents(end) frameOfNextEvent]),pk*1.05, 'rv', 'MarkerSize',3)

        xpatch = nan(4,length(onBouts));
        for j=1:length(onBouts)
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





    

end