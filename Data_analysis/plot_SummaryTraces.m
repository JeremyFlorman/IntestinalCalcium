function plot_SummaryTraces(wormdata)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
    wormdata = evalin("caller",'wormdata');
end
 
fps = 15;
nFrames = length(wormdata.bulkSignal);
bulkSignal = wormdata.bulkSignal-wormdata.backgroundSignal;
backgroundSignal  = wormdata.backgroundSignal;
autoAxialSignal = wormdata.autoAxialSignal;
area = wormdata.area;
stimTimes = wormdata.stimTimes;
velocity = wormdata.velocity;
pktraces = wormdata.peakTraces;
onFood = wormdata.onFood;
offFood = wormdata.offFood;
foodTrace = false(length(bulkSignal), 1);
if length(offFood)<length(onFood)
    offFood(end+1) = length(bulkSignal);
end

for i = 1:length(onFood)
foodTrace(onFood(i):offFood(i),1) = 1;
end



%% Plot traces
    if ~exist('time','var')
        time = linspace(0,round((nFrames)/fps/60,2),nFrames); %minutes per frame
    end
    if ~exist('pk','var')
        [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',3, 'MinPeakDistance',150);
        peakpad = fps*15;
        pktime = linspace(-15,15, peakpad*2)';
        pkmean = mean(pktraces,2,'omitnan');
    end



    figure('Position', [936 72 903 586],Color=[1 1 1])
    t = tiledlayout(4,4,'TileSpacing','compact','Padding','tight');

    % % % Bulk Signal % % %
    nexttile([1 3])
    if ~isnan(loc)
        plot(time,bulkSignal,time(loc),pk*1.01, 'rv')
    else
        plot(time,bulkSignal)
    end
    hold on
    if ~isempty(stimTimes)
        ax = gca;
        plot(time(stimTimes),ax.YLim(2)*.98,'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
    end
    ax = gca;
    plottingFoodTrace = nan(length(bulkSignal),1);
    plottingFoodTrace(foodTrace) = ax.YLim(2)*0.99;
    plot(time, plottingFoodTrace, 'Color', [0.93 0.69 0.13], 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none')
    hold off

    xlim([0 time(end)])
    %     xlabel(gca, 'Time (min)')
    ylabel(gca,'Fluorescence (a.u.)');
    title(gca, 'Whole Animal Calcium Trace')
    ax = gca;
    xt = ax.XTick;
    xtl = ax.XTickLabels;
    ax.TickLength =[0.005 0.005];
    box off

    % % % Peak Profile % % %
    nexttile([1 1])
    plot(pktime, pktraces, 'Color', [0.7 0.7 0.7])
    hold on
    plot(pktime, pkmean, 'Color', [1 0 0], 'LineWidth', 2);
    hold off
    title(gca, 'Spike Profile');
    ylabel(gca,'Fluorescence (a.u.)');
    colormap bone
    box off

    % % % Axial Signal % % %
    ax = nexttile([1 3]);
    imagesc(smoothdata(autoAxialSignal,1,'gaussian',60)'-median(backgroundSignal,'omitnan'))
    title(gca, 'Axial Calcium Trace')
    hold on
    % plot(loc,1, 'vw', 'MarkerFaceColor' ,[.4 .5 .6]);
    % if ~isempty(stimTimes)
    %     plot(stimTimes,1,'Marker', 'diamond', 'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
    % end
    ax = gca;
    plottingFoodTrace(foodTrace) = 1;
    plot(plottingFoodTrace, 'Color', [0.93 0.69 0.13], 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none')
 
    hold off
    box off

    %     xlabel('Time (min)')
    ax.XTick = xt*60*fps; %linspace(0,length(autoAxialSignal),length(xtl));
    ax.XTickLabels = xtl;
    ax.YTick = [20 size(autoAxialSignal,2)-20];
    ax.YTickLabel = {'Head', 'Tail'};
    ax.CLim =[0 50];
    colormap turbo
    ax.TickLength = [0.001 0.001];


    % % % Interval Histogram % % %
    nexttile([1 1])
    edges = 0:2:120;
    histogram(diff(loc)./fps,'BinEdges',edges);
    title(gca,'Inter-Peak Interval');
    ylim([0 10])
    xlim([0 120])
    xlabel(gca,'Time (s)');
    ylabel(gca,'Count');
    box off

    % % % Velocity  % % %
    nexttile([1 3]);
    plot(time,smoothdata(velocity,'gaussian',30))
    ax = gca;
    hold on
    if ~isempty(stimTimes)
        plot(time(stimTimes),ax.YLim(2)*.98,'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
    end
    hold off
    xlim([0 time(end)])
    title(gca, 'Velocity')
    ylabel(gca,'Steps/sec')
    %     xlabel(gca,'Time (min)')
    ax.TickLength = [0.005 0.005];
    
    box off

    % % % Peak Widths % % %
    nexttile([1 1])
    histogram(w./fps,'BinEdges', 1:15);
    ylim([0 10])
    xlim([0 15])
    title(gca,'Peak Widths');
    ylabel(gca,'Count');
    xlabel(gca,'Time (s)');
    box off
    % % % Worm Area % % %
    nexttile([1 3]);
    plot(time,smoothdata(area,'gaussian', 30))
    xlim([0 time(end)])
    title(gca, 'Worm Area')
    ylabel(gca,'Pixels')
    xlabel(gca,'Time (min)')
    ax =  gca;
    ax.TickLength = [0.005 0.005];
    box off


    % function annotateFood(data,ax)
    %      for k = 1:length(data.onFood)
    %             foodStart = data.onFood(k);
    % 
    %             if k<=length(data.offFood)
    %                 foodEnd = data.offFood(k);
    %             else
    %                 foodEnd = find(~isnan(data.bulkSignal),1,'last');
    %             end
    %             foodY = ax.YLim; 
    %             line([foodStart foodEnd], [foodY(2) foodY(2)], 'Color', [0.97 0.93 0.62], 'LineWidth', 0.5, 'LineStyle', '-', 'Marker', 'none')
    %      end
    % end

end