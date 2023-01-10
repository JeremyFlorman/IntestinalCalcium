function [mtintervalMatrix] = plotEchos(datapath, plotcontrol,labelXAxis,labelYAxis)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

intmethod = 1;  % intmethod = 1 plots the correlation between every interval. 
                % intmethod = 2 plots correlation between every other
                % interval

% peakthreshold = 750;
% plotcontrol = 1;
% datapath = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat';
% figure()

[mtdata, wtdata]= parseWormData(datapath);

mtintervalMatrix = [];
wtintervalMatrix = [];

settings = returnPlotSettings();
axlimits = settings.binedges;

peakthreshold = settings.peakthreshold;
peakwidth = settings.peakwidth;
peakdistance = settings.peakdistance;
wtcolor = settings.wtcolor; 
mtcolor = settings.mtcolor; 
mtedgecolor = settings.mtedgecolor; 
normalize = settings.normalize;
regline = 0; % add regression line? 1=yes, 0=no;

%     nidx = ~cellfun(@isempty,{mtdata.normalizedSignal});
%     mtbulksig = mtdata(nidx).normalizedSignal;


for i = 1:length(mtdata)
    if normalize == 1
        sig = mtdata(i).normalizedSignal;
    else 
        sig = mtdata(i).bulkSignal;
    end

    terpdata = fillmissing(sig, 'movmedian',100);
    [~, templocs] = findpeaks(terpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance,'MinPeakWidth',peakwidth);
    if length(templocs)>2
        mtintervals = diff(templocs)/15;
        switch intmethod
            case 1
                [mtcurrInt,mtprevInt] = getCorrelation(mtintervals);
            case 2
                [mtcurrInt, mtprevInt] = getEveryOtherCorrelation(mtintervals);
        end

        intmat = [mtcurrInt mtprevInt];
        mtintervalMatrix = vertcat(mtintervalMatrix,intmat);

    end
end


if plotcontrol ==1
%     idx = ~cellfun(@isempty,{wtdata.normalizedSignal});
%     wtbulksig = wtdata(idx).normalizedSignal;
    
    for i = 1:length(wtdata)

        if normalize == 1
            sig = wtdata(i).normalizedSignal;
        else 
            sig = wtdata(i).bulkSignal;
        end

        wtterpdata = fillmissing(sig, 'movmedian',100);
        [~, wttemplocs] = findpeaks(wtterpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance,'MinPeakWidth',peakwidth);
        if length(wttemplocs)>2
            wtintervals = diff(wttemplocs)/15;
            switch intmethod
                case 1
                    [wtcurrInt,wtprevInt] = getCorrelation(wtintervals);
                case 2
                    [wtcurrInt, wtprevInt] = getEveryOtherCorrelation(wtintervals);
            end

        end
        if ~isempty(wttemplocs)
        wtintmat = [wtcurrInt wtprevInt];
        wtintervalMatrix = vertcat(wtintervalMatrix,wtintmat);
        else
            plotcontrol = 0;
        end
    end
end

ax = gca;
ax.YLim = [axlimits(1),axlimits(end)];
ax.XLim = [axlimits(1),axlimits(end)];

if plotcontrol == 1 
    avgint = median(wtintervalMatrix,"all");
    avmint = median(mtintervalMatrix,"all");
    stdmint = std(mtintervalMatrix,1,"all");
elseif plotcontrol == 0 
    avgint = median(mtintervalMatrix,"all");
    stdint = std(mtintervalMatrix,1,"all");
end

hold on

% errorcolor = [0.1 0.1 1];
if plotcontrol == 1
    %% error bars for mutant
%     pt1 = avmint-stdmint;
%     pt2 = avmint+stdmint; 
% 
%     line([avmint avmint], [pt1 pt2],  'LineStyle', '-.','Color', errorcolor, 'LineWidth', 0.5)
%     line([pt1 pt2], [avmint avmint], 'LineStyle', '-.','Color', errorcolor, 'LineWidth', 0.5)
%         line([avmint-2 avmint+2], [pt1 pt1], 'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
%         line([avmint-2 avmint+2], [pt2 pt2], 'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
%         line([pt1 pt1], [avmint-2 avmint+2], 'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
%         line([pt2 pt2], [avmint-2 avmint+2],'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
%%

    s1 = scatter(wtintervalMatrix(:,1), wtintervalMatrix(:,2), 7,wtcolor, 'filled','MarkerEdgeColor',[0.6 0.6 0.6],...
        'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.6,'Parent',ax);
    
    s2 = scatter(mtintervalMatrix(:,1),mtintervalMatrix(:,2), 7, mtcolor, 'filled','MarkerEdgeColor',mtedgecolor,...
        'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.9,'Parent',ax);
    hold off

    if regline ==1
    ls = lsline;

    mtr = mtintervalMatrix(:,1)\mtintervalMatrix(:,2);
    wtr = wtintervalMatrix(:,1)\wtintervalMatrix(:,2);

    wid = 2;
    styl = ':';
    ls(1).Color = mtcolor;
    ls(2).Color = wtcolor;
    ls(1).LineWidth = wid;
    ls(2).LineWidth = wid;
    ls(1).LineStyle = styl;
    ls(2).LineStyle = styl;
    end

    title('Inter-interval Correlation')

lgd=legend([s1, s2], {'Control', ['\it' mtdata(1).genotype]},'Location', 'northwest');
    legend('boxoff')
%     lgd.Position = [0.7231,0.3326,0.0903,0.0454]
    if labelXAxis == 1
    xlabel('Current Interval (s)')
    end

    if labelYAxis == 1
    ylabel('Previous Interval (s)')
    end
     line(ax.XLim, [avgint avgint], 'LineStyle', ':','Color', 'k','LineWidth', 1)
    line([avgint avgint], ax.YLim, 'LineStyle', ':', 'Color', 'k','LineWidth',1) 


    lgd.String = lgd.String(1:2);

elseif plotcontrol == 0
    %% error bars for control
%     pt1 = avgint-stdint;
%     pt2 = avgint+stdint; 
% 
%     line([avgint avgint], [pt1 pt2],  'LineStyle', '-.','Color', errorcolor, 'LineWidth', 0.5)
%     line([pt1 pt2], [avgint avgint], 'LineStyle', '-.','Color', errorcolor, 'LineWidth', 0.5)
%         line([avgint-2 avgint+2], [pt1 pt1], 'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
%         line([avgint-2 avgint+2], [pt2 pt2], 'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
%         line([pt1 pt1], [avgint-2 avgint+2], 'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
%         line([pt2 pt2], [avgint-2 avgint+2],'LineStyle', '-','Color', errorcolor, 'LineWidth', 0.5)
    %%

    scatter(mtintervalMatrix(:,1),mtintervalMatrix(:,2), 7,'MarkerEdgeColor',[0.6 0.6 0.6],...
        'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.9,'Parent',ax);

    if regline == 1
    ls2 = lsline;
    wid = 2;
    styl = ':';
    ls2.Color = mtcolor;
    ls2.LineWidth = wid;
    ls2.LineStyle = styl;
    end

    title('Inter-interval Correlation')

    line(ax.XLim, [avgint avgint], 'LineStyle', ':','Color', 'k','LineWidth', 1)
    line([avgint avgint], ax.YLim, 'LineStyle', ':', 'Color', 'k','LineWidth',1)    


    if labelXAxis == 1
    xlabel('Current Interval (s)')
    end

    if labelYAxis == 1
    ylabel('Previous Interval (s)')
    end

end

end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
function [currInt, prevInt] = getCorrelation(intervals)
currInt = intervals(2:end);
prevInt = intervals(1:end-1);
end

function [currInt, prevInt] = getEveryOtherCorrelation(intervals)
intlen = [];
intidx = 1;
switch mod(length(intervals),2)
    case 0
        intlen = length(intervals);
    case 1
        intlen =length(intervals)-1;
end

for k = 2:2:intlen
    prevInt(intidx,1) = intervals(k-1);
    currInt(intidx,1) = intervals(k);
    intidx = intidx+1;
end

end
