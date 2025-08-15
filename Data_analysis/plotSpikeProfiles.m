function [] = plotSpikeProfiles(mtdata, wtdata,settings, labelXAxis,labelYAxis)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% mtdatapath = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat";
% plotcontrol = 1;
% peakthreshold = 750;
% traceylimits = [5000 12500];

if ~isfield(wtdata,'bulkSignal')
    plotcontrol = 0;
else
    plotcontrol = 1;
end


averageProfiles = 1;

traceylimits = settings.traceylimit;
secondsPrePost = settings.spikeProfileWindow;
fps = settings.framerate;
noalpha = 1;

timePreSpike = fps*secondsPrePost;
timePostSpike = fps*secondsPrePost;

% get profiles and averages per recording
if plotcontrol ==1
    wtprofiles = [wtdata(:).peakTraces];
    wtAverages = nan(timePreSpike+timePostSpike+1, length(wtdata));
    for i = 1:length(wtdata)
        if ~isempty(wtdata(i).peakTraces)
        wtAverages(:,i) = mean(wtdata(i).peakTraces, 2, "omitmissing");
        end
    end

    if averageProfiles == 1
        wtprofiles2plot = wtAverages;
    else
        wtprofiles2plot = wtprofiles;
    end
end

mtprofiles = [mtdata(:).peakTraces];
mtAverages = nan(timePreSpike+timePostSpike+1, length(mtdata));
for i = 1:length(mtdata)
    if ~isempty(mtdata(i).peakTraces)
    avg = mean(mtdata(i).peakTraces, 2, "omitmissing");
    mtAverages(1:length(avg),i) = avg;
    end
end

if averageProfiles == 1
    mtprofiles2plot = mtAverages;
else
    mtprofiles2plot = mtprofiles;
end




%% plotting
singlewidth = 0.5;
tracewidth = 1;
tracelinewidth = 1;
wtcolor = [0 0 0];
mtcolor = [0.09 0.35 0.92];

if noalpha == 1
    mttracecol = [0.7,0.85,0.99];
    wttracecol = [0.6 0.6 0.6];
else
    mttracecol = [0.7,0.85,0.99,0.6];
    wttracecol = [0.6 0.6 0.6 0.4];
end


if plotcontrol == 0
    mtcolor = wtcolor;
    mttracecol = wttracecol;

end

mtpktime = linspace(-timePreSpike/fps,timePostSpike/fps,size(mtprofiles2plot,1))';
mtpktimemat = repmat(mtpktime, [1,size(mtprofiles2plot,2)]);

hold on
if plotcontrol ==1

    wtpktime = linspace(-timePreSpike/fps,timePostSpike/fps,size(wtprofiles2plot,1))';
    wtpktimemat = repmat(wtpktime, [1,size(wtprofiles2plot,2)]);
    wtpk = max(mean(wtprofiles,2,'omitnan'));

    if ~isempty(wtpktimemat)
        plot(wtpktimemat, wtprofiles2plot, 'Color',wttracecol , 'LineWidth', singlewidth) % plot wild type individual traces
    end
end

if ~isempty(mtpktimemat)
    plot(mtpktimemat, mtprofiles2plot, 'Color',mttracecol ,'LineWidth', singlewidth); % plot mutant individual traces
    plot(mtpktime, mean(mtprofiles,2,'omitnan'), 'Color', mtcolor, 'LineWidth', tracewidth); % plot mutant trace mean
end

if plotcontrol == 1
    if~isempty(wtprofiles)
        plot(wtpktime, mean(wtprofiles,2,'omitnan'), 'Color', wtcolor, 'LineWidth', tracewidth); % plot wild type trace mean
    end
end

hold off


normalization = settings.normalize;
titles = {'Mean Fluorescence (a.u.)','\delta F/F (a.u)', 'Z-Score (s.d.)', 'Normalized Fluorescence (a.u.)'};

xlim([-timePreSpike/fps; timePostSpike/fps])
ylim(traceylimits)
title('Ca^2^+ Spike Profiles');
if labelYAxis == 1
    if strcmp(normalization, 'None') == 1
        ylabel(titles(1))
    elseif strcmp(normalization, 'Delta F/F0') == 1
        ylabel(titles(2))
    elseif strcmp(normalization, 'Z-Score') == 1
        ylabel(titles(3))
    elseif strcmp(normalization, 'Control') == 1
        ylabel(titles(4))
    end
elseif labelYAxis == 0
    yticklabels([]);
end

if labelXAxis ==1
    xlabel('Time (s)');
end

%% t-test
if plotcontrol == 1
    wtamp = max(wtprofiles);
    mtamp = max(mtprofiles);
    if ~isempty(wtamp) && ~isempty(mtamp)
        [~, p] = ttest2(mtamp,wtamp);
        if p< 0.00001
            sigval = 'p <0.00001';
        else
            sigval = ['p=' num2str(round(p,5))];
        end
        ax1 = gca();
        % text(ax1.XLim(1)+abs(ax1.XLim(1))*0.05,ax1.YLim(2)-ax1.YLim(2)*.05,sigval,'FontSize', 10)
    end
end



end


