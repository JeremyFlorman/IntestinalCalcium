function [] = plotSpikeProfiles(mtdatapath, plotcontrol,labelXAxis,labelYAxis)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% mtdatapath = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat";
% plotcontrol = 1;
% peakthreshold = 750;
% traceylimits = [5000 12500];

settings = returnPlotSettings();
traceylimits = settings.traceylimit;
peakthreshold = settings.peakthreshold;
normalize = settings.normalize;

[mtdata, wtdata] = parseWormData(mtdatapath);


wtprofiles = [];
mtpeakprofiles = [];

% midx = ~cellfun(@isempty,{mtdata.normalizedSignal});
% mtbulksig = mtdata(midx).normalizedSignal;
%

for k = 1:length(mtdata)

    if normalize == 1
        mtsig = mtdata(k).normalizedSignal;
    else
        if isfield(mtdata, 'backgroundSignal')
            background = mtdata(k).backgroundSignal;
            rawsig = mtdata(k).bulkSignal;

            mtsig = rawsig-background;
        else
            mtsig = mtdata(k).bulkSignal;
        end
    end

    terpdata = fillmissing(mtsig, 'movmedian',100);
    [pk, loc] = findpeaks(terpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',150);



    for q = 1:length(loc)
        pre = loc(q)-150;
        post = loc(q)+150;

        if pre>0 && post<= length(mtsig)
            ttrace = mtsig(pre:post);

            mtpeakprofiles = horzcat(mtpeakprofiles, ttrace);
        end
    end
end


%% Control profiles
if plotcontrol == 1

    %     widx = ~cellfun(@isempty,{wtdata.normalizedSignal});
    %     wtbulksig = wtdata(widx).normalizedSignal;

    for k = 1:length(wtdata)
        if normalize == 1
            wtsig = wtdata(k).normalizedSignal;
        else
            if isfield(wtdata, 'backgroundSignal')
                background = wtdata(k).backgroundSignal;
                rawsig = wtdata(k).bulkSignal;

                wtsig = rawsig-background;
            else
                wtsig = wtdata(k).bulkSignal;
            end
        end

        terpdata = fillmissing(wtsig, 'movmedian',100);
        [pk, loc] = findpeaks(terpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',150);


        for q = 1:length(loc)
            pre = loc(q)-150;
            post = loc(q)+150;
            if pre>0 && post<= length(wtsig)
                wttrace = wtsig(pre:post);
                wtprofiles = horzcat(wtprofiles, wttrace);
            end
        end
    end
end


%% plotting
singlewidth = 0.5;
tracewidth = 1;
tracelinewidth = 1;
wtcolor = [0 0 0];
mtcolor = [0.09 0.35 0.92];
mttracecol = [0.7,0.85,0.99,0.6];
wttracecol = [0.6 0.6 0.6 0.4];

if plotcontrol == 0
    mtcolor = wtcolor;
    mttracecol = wttracecol;
end



wtpktime = linspace(-15,15,size(wtprofiles,1))';
wtpktimemat = repmat(wtpktime, [1,size(wtprofiles,2)]);
wtpk = max(mean(wtprofiles,2,'omitnan'));

mtpktime = linspace(-15,15,size(mtpeakprofiles,1))';
mtpktimemat = repmat(mtpktime, [1,size(mtpeakprofiles,2)]);
mtpk = max(mean(mtpeakprofiles,2,'omitnan'));




hold on
if plotcontrol ==1
    if ~isempty(wtpktimemat)
        plot(wtpktimemat, wtprofiles, 'Color',wttracecol , 'LineWidth', singlewidth) % plot wild type individual traces
        line([-15; 15], [wtpk; wtpk], 'Color', wtcolor, 'LineStyle', ':',...    %line at wild type mean
            'LineWidth', tracelinewidth)
    end
end

if ~isempty(mtpktimemat)
    plot(mtpktimemat, mtpeakprofiles, 'Color',mttracecol ,'LineWidth', singlewidth); % plot mutant individual traces
    line([-15; 15], [mtpk; mtpk], 'Color', mtcolor, 'LineStyle', '--',... % line at mutant mean
        'LineWidth', tracelinewidth)

    plot(mtpktime, median(mtpeakprofiles,2,'omitnan'), 'Color', mtcolor, 'LineWidth', tracewidth); % plot mutant trace mean
end

if plotcontrol == 1
    if~isempty(wtprofiles)
        plot(wtpktime, median(wtprofiles,2,'omitnan'), 'Color', wtcolor, 'LineWidth', tracewidth); % plot wild type trace mean
    end
end

hold off




xlim([-15 15])
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
    mtamp = max(mtpeakprofiles);
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


