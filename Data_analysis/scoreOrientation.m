function  [totalScore, metrics] = scoreOrientation(prevTrace, currentTrace)
%scoreOrientation Orientation scoring for straightened worm midline signals
% PURPOSE
%   Compare frame N ("currentTrace") to frame N-1 ("prevTrace") and compute:
%       1) raw profile correlation score
%       2) spatial gradient correlation score
%       3) posterior-vs-anterior variance score
%   Then combine them into one flip score.
%
% INTERPRETATION
%   totalScore > 0   -> current frame matches previous frame better as-is
%   totalScore < 0   -> current frame matches previous frame better flipped
%
% USAGE
%   if i > 1 && ~isempty(prevTrace)
%       [totalScore, metrics] = scoreOrientation(prevTrace, currentTrace);
%       if totalScore < -0.05
%           currentTrace = fliplr(currentTrace);   % flip current row
%       end
%   end
%
% NOTES
%   - Works on row or column vectors
%   - Ignores NaNs pairwise
%   - Adds light smoothing before gradient calculation
%   - Uses fast normalized dot products instead of corr()
% -------------------------------------------------------------------------% -------------------------------------------------------------------------
% Orientation scoring for straightened worm midline signals
%
% PURPOSE
%   Compare frame N ("currentTrace") to frame N-1 ("prevTrace") and compute:
%       1) raw profile correlation score
%       2) spatial gradient correlation score
%       3) posterior-vs-anterior variance score
%   Then combine them into one flip score.
%
% INTERPRETATION
%   totalScore > 0   -> current frame matches previous frame better as-is
%   totalScore < 0   -> current frame matches previous frame better flipped
%
% USAGE INSIDE YOUR LOOP
%   if i > 1 && ~isempty(prevTrace)
%       [totalScore, metrics] = scoreOrientation(prevTrace, tt);
%       if totalScore < -0.05
%           tt = fliplr(tt);   % flip current row
%       end
%   end
%
% NOTES
%   - Works on row or column vectors
%   - Ignores NaNs pairwise
%   - Adds light smoothing before gradient calculation
%   - Uses fast normalized dot products instead of corr()
% -------------------------------------------------------------------------



% ----- force row vectors -----
prevTrace = prevTrace(:)';
currentTrace        = currentTrace(:)';

if numel(prevTrace) ~= numel(currentTrace)
    error('prevTrace and tt must have the same number of samples.');
end

% ----- parameters -----
endFrac = 0.10;     % first/last 10% for variance comparison
smoothWin = 3;      % light smoothing before gradient calculation

% ----- 1) RAW PROFILE CORRELATION SCORE -----
% Compare current frame as-is vs flipped, relative to previous frame.
sameRaw = normalizedDotPairwise(prevTrace, currentTrace);
flipRaw = normalizedDotPairwise(prevTrace, fliplr(currentTrace));

% Positive rawScore = keep as-is
% Negative rawScore = flipped version matches better
rawScore = sameRaw - flipRaw;

% ----- 2) SPATIAL GRADIENT CORRELATION SCORE -----
% Smooth a bit first so diff() is not dominated by sample-to-sample noise.
prevSm = movmean(prevTrace, smoothWin, 'omitnan');
ttSm   = movmean(currentTrace,        smoothWin, 'omitnan');

gPrev = diff(prevSm);
gSame = diff(ttSm);
gFlip = diff(fliplr(ttSm));

sameGrad = normalizedDotPairwise(gPrev, gSame);
flipGrad = normalizedDotPairwise(gPrev, gFlip);

% Positive gradScore = keep as-is
% Negative gradScore = flipped version matches better
gradScore = sameGrad - flipGrad;

% ----- 3) VARIANCE ASYMMETRY SCORE -----
% Posterior tends to be more variable/brighter in your intestinal signal.
n = numel(currentTrace);
k = max(3, round(endFrac * n));

head = currentTrace(1:k);
tail = currentTrace(end-k+1:end);

head = head(~isnan(head));
tail = tail(~isnan(tail));

if numel(head) >= 2 && numel(tail) >= 2
    varHead = var(head, 1);
    varTail = var(tail, 1);

    % Positive varScore = tail more variable than head -> likely correct
    % Negative varScore = head more variable than tail -> likely flipped
    varScore = log((varTail + eps) / (varHead + eps));
else
    varScore = 0;
end

% ----- 4) COMBINE SCORES -----
% Weight temporal continuity most strongly.
totalScore = ...
    1.0 * rawScore + ...
    1.5 * gradScore + ...
    0.5 * varScore;

% ----- output details -----
metrics.sameRaw  = sameRaw;
metrics.flipRaw  = flipRaw;
metrics.rawScore = rawScore;

metrics.sameGrad  = sameGrad;
metrics.flipGrad  = flipGrad;
metrics.gradScore = gradScore;

metrics.varScore = varScore;
metrics.totalScore = totalScore;
end


function r = normalizedDotPairwise(a, b)
% Fast Pearson-like similarity using pairwise valid samples.
% Equivalent in spirit to corr(a,b,'Rows','complete'), but faster for vectors.

a = a(:);
b = b(:);

valid = ~isnan(a) & ~isnan(b);
a = a(valid);
b = b(valid);

if numel(a) < 2
    r = 0;
    return
end

a = a - mean(a);
b = b - mean(b);

na = norm(a);
nb = norm(b);

if na < eps || nb < eps
    r = 0;
else
    r = dot(a, b) / (na * nb);
end
end

