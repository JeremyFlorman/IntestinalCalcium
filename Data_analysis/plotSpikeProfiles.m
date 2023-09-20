function [] = plotSpikeProfiles(mtdata, wtdata,settings, labelXAxis,labelYAxis)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% mtdatapath = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat";
% plotcontrol = 1;
% peakthreshold = 750;
% traceylimits = [5000 12500];

if isempty(wtdata)
    plotcontrol = 0;
else
    plotcontrol = 1;
end


traceylimits = settings.traceylimit;
secondsPrePost = settings.spikeProfileWindow;
fps = settings.framerate;
noalpha = 1;

timePreSpike = fps*secondsPrePost;
timePostSpike = fps*secondsPrePost;


if plotcontrol ==1
    wtprofiles = [wtdata(:).peakTraces];
end

mtprofiles = [mtdata(:).peakTraces];



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





mtpktime = linspace(-timePreSpike/15,timePostSpike/15,size(mtprofiles,1))';
mtpktimemat = repmat(mtpktime, [1,size(mtprofiles,2)]);
mtpk = max(mean(mtprofiles,2,'omitnan'));




hold on
if plotcontrol ==1

    wtpktime = linspace(-timePreSpike/15,timePostSpike/15,size(wtprofiles,1))';
    wtpktimemat = repmat(wtpktime, [1,size(wtprofiles,2)]);
    wtpk = max(mean(wtprofiles,2,'omitnan'));
    
    if ~isempty(wtpktimemat)
        plot(wtpktimemat, wtprofiles, 'Color',wttracecol , 'LineWidth', singlewidth) % plot wild type individual traces
        line([-timePreSpike/15; timePostSpike/15], [wtpk; wtpk], 'Color', wtcolor, 'LineStyle', ':',...    %line at wild type mean
            'LineWidth', tracelinewidth)
    end
end

if ~isempty(mtpktimemat)
    plot(mtpktimemat, mtprofiles, 'Color',mttracecol ,'LineWidth', singlewidth); % plot mutant individual traces
    line([-timePreSpike/15; timePostSpike/15], [mtpk; mtpk], 'Color', mtcolor, 'LineStyle', '--',... % line at mutant mean
        'LineWidth', tracelinewidth)

    plot(mtpktime, mean(mtprofiles,2,'omitnan'), 'Color', mtcolor, 'LineWidth', tracewidth); % plot mutant trace mean
end

if plotcontrol == 1
    if~isempty(wtprofiles)
        plot(wtpktime, mean(wtprofiles,2,'omitnan'), 'Color', wtcolor, 'LineWidth', tracewidth); % plot wild type trace mean
    end
end

hold off




xlim([-timePreSpike/fps; timePostSpike/fps])
ylim(traceylimits)
title('Ca^2^+ Spike Profiles');
if labelYAxis == 1
    ylabel('Signal Intensity (a.u.)')
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
        text(ax1.XLim(1)+abs(ax1.XLim(1))*0.05,ax1.YLim(2)-ax1.YLim(2)*.05,sigval,'FontSize', 10)
    end
end



end


