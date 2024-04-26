function [] = overlayBulkSignal(data, settings,labelXAxis)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

traceylimit = settings.traceylimit;
plotlimit = settings.tolimit;

time = 1:length(data(1).autoAxialSignal);
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
    baseline = repmat(traceylimit(1)+shift, [length(time),1]);

    %[0.7 0.2 0.4 0.7]
    hold on
    plot(time, shiftedSignal,'Color', [.9 .9 .9] ,'Marker', 'none', ...
        'LineWidth',0.75, 'LineStyle', '-')

    line(time,baseline, 'Color', [0 0 0])


end
hold off
ax = gca;
ax.YTickLabel = [];

ax.YLim = [traceylimit(1) num2plot*tracediff+traceylimit(1)];

title(['\it' data(1).genotype])
box off
ax.YColor = [1 1 1];

if labelXAxis ==1
    xlabel('Time (min)')
else
    ax.XTickLabel = [];
end

end