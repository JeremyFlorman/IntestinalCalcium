function plotEchos(mtdata,wtdata,settings,labelXAxis,labelYAxis)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

intmethod = 1;  % intmethod = 1 plots the correlation between every interval.
% intmethod = 2 plots correlation between every other
% interval

% peakthreshold = 750;
% plotcontrol = 1;
% datapath = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat';
% figure()

if isempty(wtdata)
    plotcontrol = 0; 
else
    plotcontrol = 1;
end


mtintervalMatrix = [];
wtintervalMatrix = [];


axlimits = settings.binedges;


wtcolor = settings.wtcolor;
mtcolor = settings.mtcolor;
mtedgecolor = settings.mtedgecolor;
regline = 0; % add regression line? 1=yes, 0=no;


for i = 1:length(mtdata)
mtint = mtdata(i).peakIntervals;
    if length(mtint)>2
        switch intmethod
            case 1
                [mtcurrInt,mtprevInt] = getCorrelation(mtint);
            case 2
                [mtcurrInt, mtprevInt] = getEveryOtherCorrelation(mtint);
        end

        intmat = [mtcurrInt mtprevInt];
        mtintervalMatrix = vertcat(mtintervalMatrix,intmat);

    end
end


if plotcontrol ==1
    %     idx = ~cellfun(@isempty,{wtdata.normalizedSignal});
    %     wtbulksig = wtdata(idx).normalizedSignal;

    for i = 1:length(wtdata)
        wtint = wtdata(i).peakIntervals;
        if length(wtint)>2
            switch intmethod
                case 1
                    [wtcurrInt,wtprevInt] = getCorrelation(wtint);
                case 2
                    [wtcurrInt, wtprevInt] = getEveryOtherCorrelation(wtint);
            end

        end
        if ~isempty(wtint)
            try
                wtintmat = [wtcurrInt wtprevInt];
                wtintervalMatrix = vertcat(wtintervalMatrix,wtintmat);
            catch
            end
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
    if ~isempty(wtintervalMatrix)
        s1 = scatter(wtintervalMatrix(:,1), wtintervalMatrix(:,2), 7,wtcolor, 'filled','MarkerEdgeColor',[0.6 0.6 0.6],...
            'MarkerFaceAlpha',.05,'MarkerEdgeAlpha',.6,'Parent',ax);
%             'Parent',ax);
    end

    if ~isempty(mtintervalMatrix)
        s2 = scatter(mtintervalMatrix(:,1),mtintervalMatrix(:,2), 7, mtcolor, 'filled','MarkerEdgeColor',mtedgecolor,...
            'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.9,'Parent',ax);
    end
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

    if  ~isempty(wtintervalMatrix) && ~isempty(mtintervalMatrix)

        lgd=legend([s1, s2], {'Control', ['\it' mtdata(1).genotype]},'Location', 'northwest');
        legend('boxoff')

        if labelXAxis == 1
            xlabel('Current Interval (s)')
        end

        if labelYAxis == 1
            ylabel('Previous Interval (s)')
        end

        line(ax.XLim, [avgint avgint], 'LineStyle', ':','Color', 'k','LineWidth', 1)
        line([avgint avgint], ax.YLim, 'LineStyle', ':', 'Color', 'k','LineWidth',1)
        lgd.String = lgd.String(1:2);
    end

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

    if ~isempty(mtintervalMatrix)
        scatter(mtintervalMatrix(:,1),mtintervalMatrix(:,2), 2,'MarkerEdgeColor',[0.6 0.6 0.6],...
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
