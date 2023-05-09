function [] = plotBulkSignal(datapath, plotcontrol,plotlimit)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% datapath = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat";
% plotcontrol = 1;
% peakthreshold = 750;
% traceylimit = [5000 12500];
% plotlimit = 3;


settings = returnPlotSettings();
peakthreshold = settings.peakthreshold;
peakdistance = settings.peakdistance;
peakwidth = settings.peakwidth;
traceylimit = settings.traceylimit;
normalize = settings.normalize;

[mtdata, wtdata] = parseWormData(datapath);

spikeProperties = getSpikeLocs(datapath,peakthreshold,1);
mtdata = mtdata(spikeProperties.mtsort);
wtdata = wtdata(spikeProperties.wtsort);


data = [];

switch plotcontrol
    case 0
        data = mtdata;
    case 1
        data = wtdata;
end

mpf = 1/900; 
time = mpf:mpf:length(data(1).bulkSignal)/900;
tracediff = traceylimit(2)-traceylimit(1);

if plotlimit == 0 || plotlimit>length(data)
    num2plot = length(data);
else
    num2plot = plotlimit;
end


% nidx = ~cellfun(@isempty,{data.normalizedSignal});
% bulksig = data(nidx).normalizedSignal;
% 
% 
% switch plotcontrol
%     case 0
%         bulksig = bulksig(:,mtsort);
%     case 1
%         bulksig = bulksig(:,wtsort);
% end


for i = 1:num2plot
    plotindex = num2plot-i+1; % use this to make sure the first plot is on top of axes.

    if normalize == 1
        sig = data(plotindex).normalizedSignal;
    elseif normalize == 0
        if  isfield(data, 'backgroundSignal')
            background = data(plotindex).backgroundSignal;
            rawsig = data(plotindex).bulkSignal;

            sig = rawsig-background;
        else
            sig = data(plotindex).bulkSignal;
        end
    end
    signal = fillmissing(sig, 'movmedian',100);
    [tempamp, templocs] = findpeaks(signal, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance,'MinPeakWidth',peakwidth);


    shift = tracediff*(i-1);
    shiftedSignal = signal+shift;
    shiftedamp = tempamp+shift;
    baseline = traceylimit(1)+shift;


    plot(time, shiftedSignal, 'k')
    baseLine = repmat(baseline, [length(time),1]);
    line(time,baseLine, 'Color', [0.6 0.6 0.6])
    hold on
    plot(templocs/900,shiftedamp+0.05*tracediff,'v','color' ,[0.7 0.2 0.4], 'MarkerSize',1)


end
hold off;
ax = gca;
ax.YTickLabel = [];
ax.YColor = 'none';
ax.YLim = [traceylimit(1) num2plot*tracediff+traceylimit(1)];
xlabel('Time (min)')
title(['\it' data(1).genotype])
box off


end