function warped = warp_trace_between_events(data, nPre, nMid, nPost)
%WARP_TRACE_BETWEEN_EVENTS Time-warp a kymograph between two events.
%
%   trace   : [T x P] matrix (time x position along animal)
%   idxFood : scalar index of food entry frame (1-based)
%   idxSpike: scalar index of Ca2+ wave onset frame (1-based)
%   nPre    : number of output samples before food
%   nMid    : number of output samples between food and wave
%   nPost   : number of output samples after wave
%
%   warped  : [(nPre + nMid + nPost) x P] warped trace

fps = 30;
if nargin<2
    nPre = 30*fps;
    nMid = 45*fps;
    nPost = 45*fps;
end

trace = data.autoAxialSignal;
idxFood = data.stimTimes(1);

if ~isempty(data(1).peakLoc)
    idxSpike = data.peakLoc(1);

    [T, P] = size(trace);
    origTime = (1:T)';

    % --- Pre-food segment: [1 .. idxFood] ---------------------------------
    tPreOrig = origTime(1:idxFood);
    tPreNew  = linspace(tPreOrig(1), tPreOrig(end), nPre)';
    preWarped = interp1(tPreOrig, trace(1:idxFood, :), tPreNew, 'linear', 'extrap');

    % --- Food->wave segment: [idxFood .. idxSpike] ------------------------
    tMidOrig = origTime(idxFood:idxSpike);
    tMidNew  = linspace(tMidOrig(1), tMidOrig(end), nMid)';
    midWarped = interp1(tMidOrig, trace(idxFood:idxSpike, :), tMidNew, 'linear', 'extrap');

    % --- Post-wave segment: [idxSpike .. T] -------------------------------
    tPostOrig = origTime(idxSpike:T);
    tPostNew  = linspace(tPostOrig(1), tPostOrig(end), nPost)';
    postWarped = interp1(tPostOrig, trace(idxSpike:T, :), tPostNew, 'linear', 'extrap');

    % --- Concatenate ------------------------------------------------------
    warped = [preWarped; midWarped; postWarped];
else
    disp(['Skipping Trace ' num2str(1i) ' - No Spikes Found'])
    warped = nan(nPre+nMid+nPost, 1);
end

end
