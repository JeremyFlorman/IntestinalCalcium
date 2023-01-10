function [] = overlayBulkSignal(datapath, plotcontrol,plotlimit)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% datapath = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat";
% plotcontrol = 1;
% peakthreshold = 750;
% traceylimit = [5000 12500];
% plotlimit = 3;

settings = returnPlotSettings();
peakthreshold = settings.peakthreshold; 
traceylimit = settings.traceylimit;
normalize = settings.normalize;
[mtdata, wtdata] = parseWormData(datapath);
data = [];

spikeProperties = getSpikeLocs(datapath,peakthreshold,1);
mtsort = spikeProperties.mtsort;
wtsort = spikeProperties.wtsort; 

switch plotcontrol
    case 0 
        data = mtdata(mtsort);
    case 1
        data = wtdata(wtsort);
end

% mpf = 1/900;
% time = mpf:mpf:length(data(1).bulkSignal)/900;
time = 1:9000;

tracediff = traceylimit(2)-traceylimit(1);

if plotlimit == 0 || plotlimit>length(data)
    num2plot = length(data);
else
    num2plot = plotlimit;
end

% nidx = ~cellfun(@isempty,{data.normalizedSignal});
% bulksig = data(nidx).normalizedSignal;

for i = 1:num2plot
    plotindex = num2plot-i+1; % use this to make sure the first plot is on top of axes. 


    if normalize == 1
        sig = data(plotindex).normalizedSignal;
    elseif normalize == 0
        if isfield(data, 'backgroundSignal')
            background = data(plotindex).backgroundSignal;
            rawsig = data(plotindex).bulkSignal;

            sig = rawsig-background;
        else
            sig = data(plotindex).bulkSignal;
        end
    end
    
    signal = fillmissing(sig, 'movmedian',100);


    shift = tracediff*(i-1);
    shiftedSignal = signal+shift;
    baseline = repmat(traceylimit(1)+shift, [length(time),1]);
    
    %[0.7 0.2 0.4 0.7]
    hold on
    plot(time, shiftedSignal,'Color', [.9 .9 .9 .7] ,'Marker', 'none', ...
        'LineWidth',0.75, 'LineStyle', '-')

    line(time,baseline, 'Color', [0 0 0])
    

end
hold off
ax = gca;
ax.YTickLabel = [];

ax.YLim = [traceylimit(1) num2plot*tracediff+traceylimit(1)];
xlabel('Time (min)')
title(['\it' data(1).genotype])
box off
ax.YColor = [1 1 1];



end