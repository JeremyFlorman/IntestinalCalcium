function plot_SummaryTraces(wormdatapath)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
    wormdata = evalin("caller",'wormdata');
    saveSummaryplot = 0;
else
    if isstruct(wormdatapath)
        wormdata = wormdatapath;
        saveSummaryplot = 0;
    elseif ischar(wormdatapath) || isstring(wormdatapath)
        sp = strsplit(wormdatapath, '\');
        localpath = ['C:\tmp\' sp{end}];
        copyfile(wormdatapath, localpath)
        load(localpath);
        delete(localpath);
        saveSummaryplot = 1;
    end
end

fps = 15;
patchAlpha = 1;
patchColor = [0.996 0.9400 0.7920]; %[0.93 0.69 0.13]
peakThreshold = 2.5;

nFrames = length(wormdata.bulkSignal);
bulkSignal = wormdata.bulkSignal-wormdata.backgroundSignal;
backgroundSignal  = wormdata.backgroundSignal;
autoAxialSignal = wormdata.autoAxialSignal;
area = wormdata.area;
pktraces = wormdata.peakTraces;

if isfield(wormdata, 'stimTimes')
    stimTimes = wormdata.stimTimes;
else 
    stimTimes = [];
end
if isfield(wormdata, 'velocity')
    velocity = wormdata.velocity;
end


%% Work with food patch data
trimOffFood = [];
if isfield(wormdata, 'onFood')
    onFood = wormdata.onFood;
    offFood = wormdata.offFood;
    % nFrames = length(wormdata.bulkSignal); % Need this for function
    foodTrace = false(length(bulkSignal), 1);


    if length(offFood)<length(onFood)
        offFood(end+1) = length(bulkSignal);
        trimOffFood = 1;
    end

    for i = 1:length(onFood)
        foodTrace(onFood(i):offFood(i),1) = 1;
    end

    if trimOffFood ==1
        offFood = offFood(1:end-1);
    end
end


%% Plot traces
if ~exist('time','var')
    time = linspace(0,round((nFrames)/fps/60,2),nFrames); %minutes per frame
end
if ~exist('pk','var')
    [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',peakThreshold, 'MinPeakDistance',150);
    peakpad = fps*15;
    pktime = linspace(-15,15, peakpad*2)';
    pkmean = mean(pktraces,2,'omitnan');
end

figure('Position', [936 72 903 586],Color=[1 1 1])
t = tiledlayout(4,4,'TileSpacing','compact','Padding','tight');

%% Bulk Signal % % %
nexttile([1 3])
if ~isnan(loc)
    plot(time,bulkSignal, 'k',time(loc),pk*1.05, 'rv', 'MarkerSize',3)
else
    plot(time,bulkSignal)
end
hold on
if ~isempty(stimTimes)
    ax = gca;
    plot(time(stimTimes),ax.YLim(2)*.98,'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
end
hold off

ax = gca;
xlim([0 time(end)])
ylabel(gca,'GCaMP Signal (a.u.)');
xlabel('Time (min)')
title(gca, 'Bulk Calcium Trace')
ax = gca;
xt = ax.XTick;
xtl = ax.XTickLabels;
ax.TickLength =[0.005 0.005];
box off

% Food Patches
if isfield(wormdata, 'onFood') && ~isempty(wormdata.onFood)

            boutData = computeFoodBouts(wormdata.onFood, wormdata.offFood, nFrames, [0.15 ax.YLim(2)]);
            % [patchX, patchY] = shadedFoodPatches(foodBouts, [patchYVals(i)+(tracediff*0.05), patchYVals(i+1)-(tracediff*0.05)]);
            patchXinMinutes = boutData.patchX/fps/60;
            p = patch(patchXinMinutes, boutData.patchY, patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
            uistack(p, 'bottom');

    % [patchX, patchY] = shadedFoodPatches(wormdata,ax.YLim);
    % p = patch(patchX/fps/60, patchY,  [0.93 0.69 0.13], 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
    % uistack(p, 'bottom');
end


%% Peak Profile % % %
nexttile([1 1])
plot(pktime, pktraces, 'Color', [0.7 0.7 0.7])
hold on
plot(pktime, pkmean, 'Color', [1 0 0], 'LineWidth', 2);
hold off
title(gca, 'Spike Profile');
ylabel(gca,'Fluorescence (a.u.)');
colormap bone
box off

%% Axial Signal % % %
ax = nexttile([1 3]);
imagesc(smoothdata(autoAxialSignal,1,'gaussian',60)'-median(backgroundSignal,'omitnan'))
title(gca, 'Axial Calcium Trace')
hold on
if isfield(wormdata, 'onFood')
ax = gca;
foodYVals(~foodTrace) = NaN;
foodYVals(foodTrace) = 10;
plot(foodYVals, 'Color', [0.93 0.69 0.13], 'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none')
end

hold off
box off

ax.XTick = xt*60*fps;
ax.XTickLabels = xtl;
ax.YTick = [20 size(autoAxialSignal,2)-20];
ax.YTickLabel = {'Head', 'Tail'};
ax.CLim =[0 90];
colormap turbo
ax.TickLength = [0.001 0.001];
xlabel('Time (min)')

% % % Interval Histogram % % %
nexttile([1 1])
edges = 0:2:120;
histogram(diff(loc)./fps,'BinEdges',edges);
title(gca,'Inter-Peak Interval');
% ylim([0 10])
xlim([0 120])
xlabel(gca,'Time (s)');
ylabel(gca,'Count');
box off

%% Velocity  % % %
if isfield(wormdata, 'velocity')
    nexttile([1 3]);
    plot(time,smoothdata(velocity,'gaussian',45),"Color",'k')
    ax = gca;
    hold on
    if ~isempty(stimTimes)
        plot(time(stimTimes),ax.YLim(2)*.98,'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
    end

    hold off
    xlim([0 time(end)])
    title(gca, 'Velocity')
    xlabel('Time (min)')
    ylabel(gca,'Speed (mm / sec)')
    ylim([0 0.5])
    ax.TickLength = [0.005 0.005];
    box off

    % Food Patches
    if isfield(wormdata, 'onFood')
        boutData = computeFoodBouts(wormdata.onFood, wormdata.offFood, nFrames, [0.01 ax.YLim(2)]);
        % [patchX, patchY] = shadedFoodPatches(foodBouts, [patchYVals(i)+(tracediff*0.05), patchYVals(i+1)-(tracediff*0.05)]);
        patchXinMinutes = boutData.patchX/fps/60;
        p = patch(patchXinMinutes, boutData.patchY, patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
        uistack(p, 'bottom');
    end
end

if isfield(wormdata, 'onFood') && ~isempty(wormdata.onFood)
    nexttile([1 1])

    ondur = boutData.onDur/fps/60;
    offdur = boutData.offDur/fps/60;
    binEdges = 0:2:20;
    h1 = histogram(ondur,'BinEdges',binEdges,'FaceColor',[0.93 0.69 0.13], Normalization='count');
    hold on
    h1 = histogram(offdur,'BinEdges',binEdges,'FaceColor',[0.6 0.6 0.6], Normalization='count');

 title('Time On/Off Food')
    xlabel('Bout Duration (min)')
    ylabel('Count')
    legend({'On Food', 'Off Food'})


    
else
    %% Peak Widths % % %
nexttile([1 1])
histogram(w./fps,'BinEdges', 1:15);
ylim([0 10])
xlim([0 15])
title(gca,'Peak Widths');
ylabel(gca,'Count');
xlabel(gca,'Time (s)');
box off
end

%% Worm Area % % %
nexttile([1 3]);
plot(time,smoothdata(area,'gaussian', 30))
xlim([0 time(end)])
title(gca, 'Worm Area')
ylabel(gca,'Pixels')
xlabel(gca,'Time (min)')
ax =  gca;
ax.TickLength = [0.005 0.005];
box off

% Food Patches
if isfield(wormdata, 'onFood')
    boutData = computeFoodBouts(wormdata.onFood, wormdata.offFood, nFrames, [1000 ax.YLim(2)]);
    % [patchX, patchY] = shadedFoodPatches(foodBouts, [patchYVals(i)+(tracediff*0.05), patchYVals(i+1)-(tracediff*0.05)]);
    patchXinMinutes = boutData.patchX/fps/60;
    p = patch(patchXinMinutes, boutData.patchY, patchColor, 'FaceAlpha', patchAlpha, 'EdgeColor', 'none');
    uistack(p, 'bottom');
end

if saveSummaryplot == 1
    split = strsplit(wormdatapath, '\');
    ttl = strrep(split(end), '_wormdata.mat', '');
    ttl = strrep(ttl, '_', ' ');
    ttl = strrep(ttl, '-', ' ');

    title(t,ttl);
    saveas(gcf, strrep(wormdatapath, 'wormdata.mat', 'Summary_Plots.png'))
end
end