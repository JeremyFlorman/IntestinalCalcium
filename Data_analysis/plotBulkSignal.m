function [] = plotBulkSignal(data, settings)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


traceylimit = settings.traceylimit;
fps = settings.framerate;
plotlimit = settings.tolimit;


time = linspace(0,length(data(1).bulkSignal)/fps/60,length(data(1).bulkSignal));
tracediff = traceylimit(2)-traceylimit(1);

if plotlimit == 0 || plotlimit>length(data)
    num2plot = length(data);
else
    num2plot = plotlimit;
end



for i = 1:num2plot
    plotindex = num2plot-i+1; % use this to make sure the first plot is on top of axes.
    signal = data(plotindex).bulkSignal;
    tempamp = data(plotindex).peakAmplitude;
    templocs = data(plotindex).peakLoc;

    shift = tracediff*(i-1);
    shiftedSignal = signal+shift;
    shiftedamp = tempamp+shift;
    baseline = traceylimit(1)+shift;


    plot(time, shiftedSignal, 'k')
    baseLine = repmat(baseline, [length(time),1]);
    line(time,baseLine, 'Color', [0.6 0.6 0.6])
    hold on
    plot(templocs/fps/60,shiftedamp+0.05*tracediff,'v','color' ,[0.7 0.2 0.4], 'MarkerSize',3)
    if isfield(data, 'stimTimes')
        stimTimes = data(plotindex).stimTimes;
        stimY = repmat(shift+tracediff*0.9,length(stimTimes),1);
        plot(stimTimes/fps/60,stimY,'v','color' ,[0 0 0],'MarkerFaceColor',[.8 .3 .4], 'MarkerSize',8)
    end


end
hold off;
ax = gca;
ax.YTickLabel = [];
ax.YColor = 'none';
ax.YLim = [traceylimit(1) num2plot*tracediff+traceylimit(1)];
ax.XLim = [0, max(time)];
xlabel('Time (min)')
title(['\it' data(1).genotype])
box off


end