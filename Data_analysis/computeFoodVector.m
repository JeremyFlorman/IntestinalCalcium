function [onFoodVector] = computeFoodVector(ROIs, dilationInMM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

pixelSize_mm = 0.002614186;
pxPerMm      = 1 / pixelSize_mm;        % ≈ 382.53 px / mm

dilationInPx = dilationInMM*pxPerMm;
%% Set up query points for ROI
if isfield(ROIs, 'headXY')
    hasHeadPoints = 1;
else
    hasHeadPoints = 0;
end

nPatches = numel(ROIs);

if hasHeadPoints == 1
    headXY = ROIs(1).headXY;
    midXY = ROIs(1).midXY;
    tailXY = ROIs(1).tailXY;

    nPoints = size(headXY,1);

    headIn = nan(nPoints, nPatches);
    midIn = nan(nPoints, nPatches);
    tailIn = nan(nPoints, nPatches);
else
    midXY = ROIs(1).midXY;

    nPoints = size(midXY,1);
    inPoints = nan(nPoints, nPatches);
end

% Compute points inside patches
fig = figure('Visible','off');
for i =1:nPatches
    c = drawcircle('Center', ROIs(i).Center,'Radius', ...
        ROIs(i).Radius+dilationInPx, 'Color',[0.9059    0.1608    0.5412], ...
        FaceAlpha=0, LineWidth=1.5, MarkerSize=1); % dilate circle by 1mm

    if hasHeadPoints == 1
        headIn(1:nPoints, i) = inROI(c, headXY(:,1), headXY(:,2));
        midIn(1:nPoints, i) = inROI(c, midXY(:,1), midXY(:,2));
        tailIn(1:nPoints, i) = inROI(c, tailXY(:,1), tailXY(:,2));
    else
        inPoints(1:nPoints, i) = inROI(c, midXY(:,1), midXY(:,2));
    end
end

close(fig)

if hasHeadPoints == 1
    onFoodVector = any(headIn,2) | any(midIn,2) | any(tailIn,2);
else
    onFoodVector = any(inPoints, 2);
end

end