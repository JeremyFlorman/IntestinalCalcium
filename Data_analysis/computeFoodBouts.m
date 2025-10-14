function [boutData] = computeFoodBouts(onFood, offFood, totalFrames, yRange)
% computeFoodBouts_withPatches
%   Pair each on/off event with the next opposite event after it (or totalFrames).
%   Also returns patch polygons (x,y) for each on-food bout for plotting.
%
% Inputs:
%   onFood, offFood : vectors (frame indices) of enters/leaves (chronological)
%   totalFrames     : scalar, total length of recording (assumed known)
%   yRange (opt)    : [ymin ymax] for patches; defaults to [0 1]
%
% Outputs:
%   onBouts  Nx2 [start end]
%   offBouts Mx2 [start end]
%   onDur, offDur : durations in frames
%   onPatches : struct array with fields .x .y for patch coords (same order as onBouts)

if nargin < 4 || isempty(yRange)
    yRange = [0 1];
end

% make column vectors and ensure sorted & unique (defensive)
onFood  = unique(sort(onFood(:)));
offFood = unique(sort(offFood(:)));

% --- build onBouts by pairing each onStart with first off > onStart
if isempty(onFood)
    onBouts = zeros(0,2);
    onDur   = zeros(0,1);
    onPatches = struct('x',{},'y',{});
else
    nOn = numel(onFood);
    onStart = onFood;
    onEnd = zeros(nOn,1);
    for i = 1:nOn
        idx = find(offFood > onStart(i), 1, 'first');
        if ~isempty(idx)
            onEnd(i) = offFood(idx);
        else
            onEnd(i) = totalFrames;   % no off after this on => runs to end
        end
    end
    % remove any zero-length or reversed bouts (defensive)
    valid = onEnd > onStart;
    onBouts = [onStart(valid) onEnd(valid)];
    onDur = onBouts(:,2) - onBouts(:,1);
end

% --- build offBouts by pairing each offStart with first on > offStart
% include initial off bout that starts at 1
offStartCandidates = [1; offFood];
nOff = numel(offStartCandidates);
offEnd = zeros(nOff,1);
for i = 1:nOff
    s = offStartCandidates(i);
    idx = find(onFood > s, 1, 'first');
    if ~isempty(idx)
        offEnd(i) = onFood(idx);
    else
        offEnd(i) = totalFrames;
    end
end
validOff = offEnd > offStartCandidates;
offBouts = [offStartCandidates(validOff) offEnd(validOff)];
offDur = offBouts(:,2) - offBouts(:,1);

boutData.onBouts = onBouts;
boutData.offBouts = offBouts;
boutData.onDur = onDur;
boutData.offDur = offDur; 

for k = 1:size(onBouts,1)
    s = onBouts(k,1); e = onBouts(k,2);
    boutData.patchX(1:4, k) = [s e e s];
    boutData.patchY(1:4, k) = [yRange(1) yRange(1) yRange(2) yRange(2)];
end
end