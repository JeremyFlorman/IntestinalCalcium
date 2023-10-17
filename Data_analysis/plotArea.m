function [] = plotArea(data, settings)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%this is a test
traceylimit = settings.traceylimit;
fps = settings.framerate;
plotlimit = settings.tolimit;


time = linspace(0,length(data(1).bulkSignal)/fps/60,length(data(1).bulkSignal));
tracediff = 1000;

if plotlimit == 0 || plotlimit>length(data)
    num2plot = length(data);
else
    num2plot = plotlimit;
end


for i = 1:num2plot
    plotindex = num2plot-i+1; % use this to make sure the first plot is on top of axes.
    
    area = rescale(smoothdata(data(plotindex).area, 'gaussian', 30));
    templocs = data(plotindex).peakLoc;

    
    shiftedSignal = area+(i-1);

    plot(time, shiftedSignal, 'k')

    baseLine = repmat(i-1, [length(time),1]);
    line(time,baseLine, 'Color', [0.6 0.6 0.6])

    locYs = repmat(i-.06,length(templocs),1);

    hold on
    plot(templocs/fps/60,locYs,'v','color' ,[0.7 0.2 0.4], 'MarkerSize',5)
    if isfield(data, 'stimTimes')
        stimTimes = data(plotindex).stimTimes;
        stimY = repmat(i,length(stimTimes),1);
        plot(stimTimes/fps/60,stimY,'v','color' ,[0 0 0],'MarkerFaceColor',[.8 .3 .4], 'MarkerSize',8)
    end


end
hold off;
ax = gca;
ax.YTickLabel = [];
ax.YColor = 'none';
ax.YLim = [0 length(data)];
ax.XLim = [0, max(time)];
xlabel('Time (min)')
title(['\it' data(1).genotype])
box off


end