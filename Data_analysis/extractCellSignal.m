function [int1Signal, int9Signal] = extractCellSignal(inputData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isstruct(inputData)
    axialSignal = inputData(1).autoAxialSignal;
elseif ismatrix(inputData)
    axialSignal = inputData;
elseif iscell(inputData)
    axialSignal = inputData{:};
end

%% Int1 and Int9 signal
dataSize = size(axialSignal);
segmentSize = round(dataSize(2)*0.025);
lengthChunk = floor(dataSize(2)/3);


% Int1 start and end points
anteriorMean = mean(axialSignal(:,1:lengthChunk),1,'omitmissing');
[~,antLoc] = findpeaks(rescale(anteriorMean), NPeaks=1, SortStr="descend");

antStart = antLoc-segmentSize;
antEnd = antLoc+segmentSize;
if antStart<1
    antStart = 1;
end

% Int9 start and end points
posteriorMean = mean(axialSignal(:,end-lengthChunk:end),1,'omitmissing');
[~,postLoc] = findpeaks(rescale(posteriorMean), NPeaks=1, SortStr="descend");

postLoc = postLoc+(lengthChunk*2); % account for the anterior two thirds
postStart = postLoc-segmentSize;
postEnd = postLoc+segmentSize;
if postEnd > dataSize(2)
    postEnd = dataSize(2);
end

    %% plot int1/int9 identification
    % plot(mean(inputData(i).autoAxialSignal,1,'omitmissing'))
    % ax = gca;
    % lims = ax.YLim;
    % line([antLoc-segmentSize antLoc-segmentSize], [lims(1) lims(2)])
    % line([antLoc+segmentSize antLoc+segmentSize], [lims(1) lims(2)])
    % line([postLoc-segmentSize postLoc-segmentSize],  [lims(1) lims(2)])
    % line([postLoc+segmentSize postLoc+segmentSize],  [lims(1) lims(2)])
    %

int1Signal = mean(axialSignal(:,antStart:antEnd),2, 'omitmissing');
int9Signal = mean(axialSignal(:,postStart:postEnd),2, 'omitmissing');
end