function [int1Signal, int9Signal] = extractCellSignal(inputData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Int1 and Int9 signal
dataSize = size(inputData.autoAxialSignal);
segmentSize = round(dataSize(2)*0.05);
lengthChunk = floor(dataSize(2)/3);


% Int1 start and end points
anteriorMean = mean(inputData.autoAxialSignal(:,1:lengthChunk),1,'omitmissing');
[~,antLoc] = findpeaks(rescale(anteriorMean), NPeaks=1, SortStr="descend");

antStart = antLoc-segmentSize;
antEnd = antLoc+segmentSize;
if antStart<1
    antStart = 1;
end

% Int9 start and end points
posteriorMean = mean(inputData.autoAxialSignal(:,end-lengthChunk:end),1,'omitmissing');
[~,postLoc] = findpeaks(rescale(posteriorMean), NPeaks=1, SortStr="descend");

postLoc = postLoc+(lengthChunk*2); % account for the anterior two thirds
postStart = postLoc-segmentSize;
postEnd = postLoc+segmentSize;
if postEnd > dataSize(2)
    postEnd = dataSize(2);
end

%% plot int1/int9 identification
% plot(mean(inputData.autoAxialSignal,1,'omitmissing'))
% ax = gca;
% lims = ax.YLim;
% line([antLoc-10 antLoc-10], [lims(1) lims(2)])
% line([antLoc+10 antLoc+10], [lims(1) lims(2)])
% line([postLoc-10 postLoc-10],  [lims(1) lims(2)])
% line([postLoc+10 postLoc+10],  [lims(1) lims(2)])
%%

int1Signal = mean(inputData.autoAxialSignal(:,antStart:antEnd),2, 'omitmissing');
int9Signal = mean(inputData.autoAxialSignal(:,postStart:postEnd),2, 'omitmissing');
end