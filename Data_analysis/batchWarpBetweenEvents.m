data = wildtypenoFood0minData;

fps = 30;

nPre = 30*fps;
nMid = 45*fps;
nPost = 45*fps;

int1Matrix = nan(nPre+nMid+nPost, numel(data));
int9Matrix = nan(nPre+nMid+nPost, numel(data));

for i = 1:length(data)
    warped = warp_trace_between_events(data(i));
    
    if ~all(isnan(warped))
        [int1, int9] = extractCellSignal(warped);
        int1Matrix(:, i) = int1;
        int9Matrix(:, i) = int9;
    end

end