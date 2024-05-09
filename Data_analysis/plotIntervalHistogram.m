function [] = plotIntervalHistogram(mtdata, wtdata,settings, showlegend, labelXAxis,labelYAxis)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% mtdatapath = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat";
% plotcontrol = 1;
% peakthreshold = 750;
% mtcolor = [0.09 0.35 0.92];


loclinewidth = 1.5;
mtcolor = settings.mtcolor;
binedges = settings.binedges;
wtcolor = [0 0 0];
noalpha = 1;

if isempty(wtdata)
    plotcontrol = 0;
else 
    plotcontrol = 1;
end


mtint = vertcat(mtdata(:).peakIntervals);



if plotcontrol == 1
    wtint = vertcat(wtdata(:).peakIntervals);

    if noalpha ==1
        h1=histogram(wtint,binedges,'FaceColor', [0.7 0.7 0.7],'FaceAlpha',...
            1, 'EdgeAlpha', 1, 'Normalization','probability');
        hold on

        h2 =  histogram(mtint,binedges,'FaceColor', mtcolor,'FaceAlpha',...
            1, 'EdgeAlpha', 1, 'Normalization','probability');
    else
        h1=histogram(wtint,binedges,'FaceColor', [0.7 0.7 0.7],'FaceAlpha',...
            0.4, 'EdgeAlpha', 0.4, 'Normalization','probability');
        hold on

        h2 =  histogram(mtint,binedges,'FaceColor', mtcolor,'Normalization','probability');
    end
    ax1 = gca;
    ax1.YLim = [0 0.5];


    line([mean(mtint); mean(mtint)], [ax1.YLim(1) ax1.YLim(2)*.8], 'Color', mtcolor, 'LineStyle', '--',...
        'LineWidth', loclinewidth)
    line([mean(wtint); mean(wtint)], [ax1.YLim(1) ax1.YLim(2)*.8], 'Color', wtcolor, 'LineStyle', ':',...
        'LineWidth', loclinewidth)

    hold off

    if labelYAxis == 1
        ylabel('Probability');
    end

    if labelXAxis == 1
        xlabel('Interval (s)')
    end


    title('Inter-spike Interval')

    if showlegend == 1
        try
            legend([h1, h2], {'Control', ['\it' mtdata(1).genotype]});
            legend('boxoff')
        catch
        end
    end
    % if ~isempty(mtint) && ~isempty(wtint)
    %     [~, p] = ttest2(mtint,wtint);
    %     if p< 0.00001
    %         sigval = 'p <0.00001';
    %     else
    %         sigval = ['p=' num2str(round(p,5))];
    %     end
    %     text(0,ax1.YLim(2)*0.75,sigval,'FontSize', 10)
    % end
elseif plotcontrol == 0

    histogram(mtint,binedges,'FaceColor', [0.7 0.7 0.7],'FaceAlpha',0.4,...
        'Normalization','probability');
    hold on
    ax2 = gca;
    ax2.YLim = [0 0.5];
    line([mean(mtint); mean(mtint)], [ax2.YLim(1) ax2.YLim(2)*.8],...
        'Color', wtcolor, 'LineStyle', '--','LineWidth', loclinewidth)

    hold off

    if labelYAxis == 1
        ylabel('Probability');
    end

    if labelXAxis == 1
        xlabel('Interval (s)')
    end

    title('Inter-spike Interval')

    if showlegend == 1
        legend(['\it' mtdata(1).genotype])
        legend('boxoff')
    end

end
box off
ax = gca;
ax.XLim = [binedges(1) binedges(end)];
end