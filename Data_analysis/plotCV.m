function [] = plotCV(mtdata,wtdata, settings,labelYAxis)
%UNTITLED2 Summary of this function goes here
%   optional inputs are 1: threshold used for detecting spikes. 2: whether
%   to plot the CV within this function or just output the values.

% mtdatapath = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat';
% peakthreshold = 750;
% plotcontrol = 0;



if isempty(wtdata)
    plotcontrol = 0;
else
    plotcontrol = 1;
end


mtcolor = settings.mtcolor;
mtedgecolor = settings.mtedgecolor;


wtcv = nan(length(wtdata),1);
mtcv = nan(length(mtdata),1);


for i = 1:length(mtdata)
tempints = mtdata(i).peakIntervals;
    if length(tempints)>2
        mtcv(i) = std(tempints)/mean(tempints,'omitnan')*100;
    end
end


if plotcontrol == 1
    for i = 1:length(wtdata)
        tempints = wtdata(i).peakIntervals;
        if length(tempints)>2
            wtcv(i) = std(tempints)/mean(tempints,'omitnan')*100;
        end
    end
end


if plotcontrol == 0
    boxplot(mtcv)
    hold on
    sz = repmat(10,length(mtcv),1);
    scatter(ones(length(mtcv),1),mtcv,sz,'MarkerEdgeColor',[0 0 0])
    hold off
    ylim([0 100])
    ax = gca;
    ax.TickLabelInterpreter = 'tex';
    ax.XTickLabel = ['\it' mtdata(1).genotype];
    title('Coefficient of Variation')
    if labelYAxis == 1
        ylabel('% of Mean')
    end
    box off

elseif plotcontrol == 1
    cvmat = NaN(max(length(wtcv),length(mtcv)),2);
    cvmat(1:length(wtcv),1) = wtcv;
    cvmat(1:length(mtcv),2) = mtcv;

    boxplot(cvmat)
    hold on
    scatter(ones(length(cvmat),1),cvmat(:,1),10,'MarkerEdgeColor',[0 0 0])
    scatter(repmat(2,length(cvmat),1),cvmat(:,2), 10, mtcolor, 'filled','MarkerEdgeColor',mtedgecolor)
    hold off


    ylim([0 100])
    ax = gca;
    ax.TickLabelInterpreter = 'tex';
    ax.XTickLabel = {'Control', 'Mutant'};% ['\it' mtdata(1).genotype]};
    title('Coefficient of Variation')
    if labelYAxis == 1
        ylabel('% of Mean')
    end
    box off




    [~, p] = ttest2(mtcv,wtcv);
    if p< 0.00001
        sigval = 'p <0.00001';
    else
        sigval = ['p=' num2str(round(p,5))];
    end

    text(ax.XLim(1)+abs(ax.XLim(1))*0.075,ax.YLim(2)-ax.YLim(2)*.05,sigval,'FontSize', 10)


end



end