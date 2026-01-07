function [] = plotBulkSignal(data, settings, labelXAxis)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%this is a test
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


patchYVals = traceylimit(1):traceylimit(2):traceylimit(2)*num2plot+1;
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
        stimY = repmat(shift+tracediff*0.75,length(stimTimes),1);
        plot(stimTimes/fps/60,stimY,'v','color' ,[0 0 0],'MarkerFaceColor',[.8 .3 .4], 'MarkerSize',5)
    end

    %% Annotate pBocs
    if isfield(data, 'pBoc')
        bocTimes = data(plotindex).pBoc;
        bocY = shift+tracediff*0.75;
        for k = 1:length(bocTimes)
            text(bocTimes(k)/fps/60, bocY, 'P', 'Color',[0 0 0], 'FontSize', 5)
        end
    end




    
    if settings.annotateFood == 1
        if isfield(data(plotindex), 'onFood') && ~isempty(data(plotindex).onFood)
            yRange =  [patchYVals(i)+(tracediff*0.05), patchYVals(i+1)-(tracediff*0.05)];
            numFrames = find(~isnan(signal), 1, "last");

            boutData = computeFoodBouts(data(plotindex).onFood, data(plotindex).offFood, numFrames, yRange);
            % [patchX, patchY] = shadedFoodPatches(foodBouts, [patchYVals(i)+(tracediff*0.05), patchYVals(i+1)-(tracediff*0.05)]);
            patchXinMinutes = boutData.patchX/fps/60;
            p = patch(patchXinMinutes, boutData.patchY,  [0.93 0.69 0.13], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
            uistack(p, 'bottom');

            % foodStart = data(plotindex).onFood;
            % if isempty(foodStart)
            %     foodStart = 1;
            % end
            % foodEnd  = nan(length(foodStart),1);
            % for k = 1:length(foodStart)
            %     if k<=length(data(plotindex).offFood)
            %         foodEnd(k) = data(plotindex).offFood(k);
            %     else
            %         foodEnd(k) = find(~isnan(data(plotindex).bulkSignal),1,'last');
            %     end
            % end
            % foodY = shift+tracediff*0.9;
            % foodX = [foodStart/fps/60 foodEnd/fps/60];
            % for j = 1:size(foodX,1)
            %     line(foodX(j,:), [foodY foodY], 'Color', [0.93 0.69 0.13], 'LineWidth', 0.5, 'LineStyle', '-', 'Marker', 'o')
            % end
        end
    end



end
hold off;
ax = gca;
ax.YTickLabel = [];
ax.YColor = 'none';
ax.YLim = [traceylimit(1) num2plot*tracediff+traceylimit(1)];
ax.XLim = [0, max(time)];
title(['\it' data(1).genotype])
box off

if labelXAxis ==1
    xlabel('Time (min)')
else
    ax.XTickLabel = [];
end

end