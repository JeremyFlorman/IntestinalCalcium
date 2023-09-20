function [] = genotypeSummaryPlot(dataPath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Settings
axylimits = [3000, 40000];
traceylimits = [5000 12500];
peakthreshold =750;
maxNumTraces = 0; % # of traces to plot. set to 0 to plot all traces
controlname = 'wildtype'; % name of control dataset
                            %(e.g. wildtype_mergedData.mat would be 'wildtype')
%%

if nargin <1
               datapath= uigetdir('C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\')
%     datapath= 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\egl-19(gf)';
end

datadir = dir([datapath '\**\*mergedData.mat']);
[genotypes, wtidx] = getGenotypes({datadir(:).name}, controlname); %'EmptyVec'
wt = load(fullfile(datadir(wtidx).folder, datadir(wtidx).name));
wtdata = wt.wormdata;
axmean = {};


for j =  1:length(genotypes)
    wd =  load(fullfile(datadir(j).folder, datadir(j).name));
    wormdata = wd.wormdata;
    savename = [datapath '\' genotypes{j} '_DataPlots.png']
    BSsavename = [datapath '\' genotypes{j} '_BulkSignal.png'];
    ASsavename = [datapath '\' genotypes{j} '_AxialSignal.png'];
    Echosavename = [datapath '\'  genotypes{j} '_Interval Correlation.png'];
    axsavename = [datapath '\AxialMean.png'];
    CIsavename = [datapath '\CI.png'];
    PRsavename = [datapath '\PropagationRate.png'];
    GTsavename = [datapath '\BulkTraces.png'];
    ISIsavename = [datapath '\InterSpikeInterval.png'];
    
    
    
    if j == wtidx
        disp('this is the control data');
        iswildtype = 1;
        [axialMean, CI, propagationRate, axialWindow, mtprofiles] = plotTraces(wormdata,wtdata,iswildtype, genotypes{j},...
            axylimits, traceylimits, peakthreshold, maxNumTraces);
    else
        iswildtype = 0;
        [axialMean, CI, propagationRate,axialWindow, mtprofiles] = plotTraces(wormdata,wtdata,iswildtype,genotypes{j},...
            axylimits, traceylimits, peakthreshold,maxNumTraces);
        
        disp('Look at these mutants');
    end
    exportgraphics(gcf, savename, 'Resolution', 300);
    axmean(j) = {axialMean};
    groupCI(j) = {CI};
    groupPR(j) = {propagationRate};
    groupTraces(j) = {mtprofiles};
    
    figure('Position', [72.2000 171.4000 330.4000 542.4000])
    plotBulkSignal2(wormdata)
    exportgraphics(gcf, BSsavename, 'Resolution', 300);

    figure('Position', [488 70.6000 465 691.4000]);
    plotAxialSignal(wormdata)
    exportgraphics(gcf, ASsavename, 'Resolution', 300);

    figure();
    plotEchos(wormdata,peakthreshold)
    exportgraphics(gcf, Echosavename, 'Resolution', 300);
end
axylimits = [3000, 50000];


plotGroupAxialMean(axmean,axsavename,genotypes,axylimits,axialWindow);
% plotGroupCI(groupCI, CIsavename, genotypes);
% plotGroupPR(groupPR, PRsavename, genotypes, axialWindow);
plotGroupTraces(groupTraces, GTsavename, genotypes, traceylimits, wtidx)
plotGroupISI(groupCI, ISIsavename, genotypes, wtidx);


close all
end


function [gtypes, wtidx] = getGenotypes(data, controlname)
gtypes = cell(length(data), 1);
wtidx = [];
for i = 1:length(data)
    r = regexp(data{i}, '_', 'split');
    r = r{1};
    gtypes(i) = {r};
    if strcmp(r, controlname) == 1
        wtidx = i;
    end
    
end
end
function [axialMean, CI, propagationRate,axialWindow, mtprofiles] = plotTraces(wormdata, wtdata, iswildtype, genotype,ylimits,traceylimits,peakthreshold,maxNumTraces)
if maxNumTraces == 0
    maxNumTraces = length(wormdata);
end

hF = figure('Position', [38.6000 93 1.4568e+03 660.8000]);
t = tiledlayout(maxNumTraces*2, 3, 'TileSpacing','Compact','Padding','Compact');
CI = [];
sumPeak = [];
mtpeakprofile = [];
mtloc = [];
wtloc = [];
wtpeakprofile = [];
axialSpikes = {};
resamplelen = 150;
CIylimits = [1*10^6, 4*10^6];





%%%%%%%% Re-extract spike Locations and Profiles/Traces for WT %%%%%%%%%%%%

for wi = 1:length(wtdata)
    
    wtterpdata = fillmissing(wtdata(wi).bulkSignal, 'movmedian',100);
    [wtpk, wloc] = findpeaks(wtterpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',150);
    
    
    %%%%% Re-Extract Spike profiles %%%%%%%%%%
    for wq = 1:length(wloc)
        pre = wloc(wq)-150;
        post = wloc(wq)+150;
        if pre>0 && post<= length(wtdata(wi).bulkSignal)
            wttrace = wtdata(wi).bulkSignal(pre:post);
            if isempty(wtpeakprofile)
                wtpeakprofile = wttrace;
            elseif length(wtpeakprofile) == length(wttrace)
                wtpeakprofile = horzcat(wtpeakprofile, wttrace);
            end
        end
    end
    
    
    if isempty(wtloc) && ~isempty(wloc)
        wtloc = diff(wloc)/15;
    elseif ~isempty(wloc)
        wtloc = vertcat(wtloc, diff(wloc)/15);
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

for k = 1:length(wormdata)
    %% plot traces
    terpdata = fillmissing(wormdata(k).bulkSignal, 'movmedian',100);
    [pk, loc] = findpeaks(terpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',150);
    
    sumpeak = NaN(length(loc), 1);
    
    %%%%% Re-Extract Spike profiles and axial signal for mutant %%%%%%%%%%
    for q = 1:length(loc)
        pre = loc(q)-150;
        post = loc(q)+150;
        
        if pre>0 && post<= length(wormdata(k).bulkSignal)
            ttrace = wormdata(k).bulkSignal(pre:post);
            sumpeak(q,1) = sum(wormdata(k).bulkSignal(pre:post),'omitnan'); % for constipation index
            if isempty(mtpeakprofile)
                mtpeakprofile = ttrace;
            elseif length(mtpeakprofile) == length(ttrace)
                mtpeakprofile = horzcat(mtpeakprofile, ttrace);
            end
            
            %% get axial spikes and resample
            
            taxpre = loc(q)-75;
            taxpost = loc(q)+75;
            
            
            if taxpre>0 && taxpost<= length(wormdata(k).autoAxialSignal)
                tax = wormdata(k).autoAxialSignal(taxpre:taxpost,:);
                for rs = 1:size(tax,1)
                    rsamp = resample(tax(rs,:), resamplelen,length(tax(rs,:)));
                    rtax(rs,1:resamplelen) = rsamp;
                end
                
                if isempty(axialSpikes) && ~isempty(rtax)
                    axialSpikes(1) = {rtax};
                elseif ~isempty(rtax)
                    axialSpikes(length(axialSpikes)+1) = {rtax};
                end
            end
        end
        
        
    end
    
    
    
    %%%%%%%%%%%%% Re-extract spike Locations for mutant %%%%%%%%%%%%%%%%%%%%%%%
    if isempty(mtloc) && ~isempty(loc)
        mtloc = diff(loc)/15;
    else
        mtloc = vertcat(mtloc, diff(loc)/15);
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
    
    if length(pk) >2
        tempint = diff(loc)/15;
        %         temppk = pk(2:end)/mean(pk, 'omitnan');
        temppk = sumpeak(2:end,1);
        if isempty(CI)
            CI = [tempint temppk];
        else
            CI = vertcat(CI, [tempint temppk]);
        end
    end
    
    %     mtterpdata = fillmissing(wormdata(k).bulkSignal, 'movmedian',100);
    %     [mtpk, mtloc] = findpeaks(mtterpdata, 'MinPeakProminence', threshold, 'MinPeakDistance',150);
    
    time = linspace(0,10,9000);
    locinmin = loc/(15*60);
    wtcolor = [0 0 0];
    mtcolor = [0.09 0.35 0.92];
    
    if iswildtype == 1
        bulkcolor = wtcolor;
    elseif iswildtype == 0
        bulkcolor = mtcolor;
    end
    
    
    if k<=maxNumTraces
        nexttile([1,3])
        plot(time,fillmissing(wormdata(k).bulkSignal, 'movmedian',10),...
            'Color', bulkcolor, 'LineWidth', 1.5)
        xlim([0 10])
        if ~isempty(pk)
            hold on
            plot(locinmin, pk*1.15, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.2 0.2 0.2]);
            hold off
        end
        ylim(traceylimits)
        shutOffLabels(gca);
        
        %% axial signal
        nexttile([1,3])
        imagesc(wormdata(k).autoAxialSignal', ylimits) %'Interpolation',"bilinear", 'Border', 'tight')
        colormap parula
        shutOffLabels(gca);
        drawnow
    end
end
%frame = frame2im(getframe(gcf));
frame = print('-RGBImage', '-r300');
close(gcf)


%% compute mean axial spike  and Propagation rate %%%%%%%%%%%%%%%%%%%%%%%%%
axialMean = NaN(size(axialSpikes{1},1),size(axialSpikes{1},2), length(axialSpikes));
propagationRate = NaN(length(axialSpikes),1);
% figure();
for b = 1:length(axialSpikes)
    axialMean(:,:, b) = axialSpikes{b};
    thisSpike = axialSpikes{b};
    chunk = floor(size(thisSpike,2)/4);
    
    headSpike = mean(thisSpike(:,1:chunk),2, 'omitnan');
    tailSpike = mean(thisSpike(:,end-(chunk-1):end),2,'omitnan');
%     [hmax,hidx] = max(headSpike,[],'omitnan');
%     [tmax,tidx] = max(tailSpike, [], 'omitnan');
%     plot(1:length(headSpike), headSpike, 1:length(tailSpike), tailSpike);
%     line([hidx hidx], [0 hmax+100])
%     line([tidx tidx], [0 tmax+100])
    try
        tpts = findchangepts(tailSpike);
        hpts = findchangepts(headSpike);
        propagationRate(b) = 1/((hpts-tpts)/15);
    catch
    end
%     drawnow()
%     pause(0.01)
end

% histogram(propagationRate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


axialMean =  mean(axialMean, 3, 'omitnan');
%%
figure('Position', [38.6000 65.8000 1.2464e+03 688])
t2 = tiledlayout(4,5);
nexttile([3 4])
imshow(frame(40:end-40, 100:end-85,:), 'Border', 'tight')
shutOffLabels(gca)
xlabel('Time (min)');
ylabel('Signal Intensity (a.u.)');

%% plot mean axial signal
win = ((taxpost-taxpre)/15)/2;
axialWindow = win;
nexttile([4 1])
imagesc(axialMean',ylimits)
colormap parula
title('Mean Axial Spike');
ca = gca;


ca.XTick = linspace(1,size(axialMean,1), 7);
ca.XTickLabels = floor(linspace(-win,win,7));
xlabel('Time (s)')

ca.YTick = [ca.YTick(1) ca.YTick(end)];
ca.YTickLabels = {'Head', 'Tail'};
ca.FontWeight = 'bold';



nexttile([1 1])
wtprofiles =  wtpeakprofile; %cell2mat({wtdata(:).peakTraces});
wtpktime = linspace(-15,15,length(wtprofiles))';
wtpktimemat = repmat(wtpktime, [1,size(wtprofiles,2)]);
wtpk = max(mean(wtprofiles,2,'omitnan'));



mtprofiles = mtpeakprofile; %cell2mat({wormdata(:).peakTraces});
mtpktime = linspace(-15,15,length(mtprofiles))';
mtpktimemat = repmat(mtpktime, [1,size(mtprofiles,2)]);
mtpk = max(mean(mtprofiles,2,'omitnan'));


singlewidth = 0.5;
tracewidth = 2;
tracelinewidth = 1.5;
loclinewidth = 1.5;

hold on
plot(wtpktimemat, wtprofiles, 'Color', [0.6 0.6 0.6], 'LineWidth', singlewidth) % plot wild type individual traces
line([-15; 15], [wtpk; wtpk], 'Color', wtcolor, 'LineStyle', ':',...    %line at wild type mean
    'LineWidth', tracelinewidth)

if iswildtype==0
    plot(mtpktimemat, mtprofiles, 'Color', [0.7,0.85,0.99],'LineWidth', singlewidth); % plot mutant individual traces
    line([-15; 15], [mtpk; mtpk], 'Color', mtcolor, 'LineStyle', '--',... % line at mutant mean
        'LineWidth', tracelinewidth)
end

plot(wtpktime, mean(wtprofiles,2,'omitnan'), 'Color', wtcolor, 'LineWidth', tracewidth); % plot wild type trace mean
if iswildtype ==0
    plot(mtpktime, mean(mtprofiles,2,'omitnan'), 'Color', mtcolor, 'LineWidth', tracewidth); % plot mutant trace mean
end


xlim([-15 15])
ylim(traceylimits)
title('Ca^2^+ Spike Profiles');
ylabel('Signal Intensity (a.u.)')
xlabel('Time (s)');
hold off

%% inter spike interval
nexttile([1 1])

wtlocs = [];
mtlocs = [];
for w = 1:length(wtdata)
    tl = diff(wtdata(w).peakLoc)/15;
    wtlocs = vertcat(wtlocs, tl);
end

for m = 1:length(wormdata)
    mtl = diff(wormdata(m).peakLoc)/15;
    mtlocs = vertcat(mtlocs, mtl);
end



wtlocs = wtloc;
mtlocs = mtloc;

if iswildtype == 1
    histogram(wtlocs,0:5:250,'FaceColor', wtcolor);
    %     ylim([0 25])
    %     xlim([0 250])
    line([median(wtlocs); median(wtlocs)], [0; 25], 'Color', wtcolor, 'LineStyle', ':',...
        'LineWidth', loclinewidth)
    
    
elseif iswildtype == 0
    hold on
    h1=histogram(wtlocs,0:5:250,'FaceColor', [0.7 0.7 0.7],'FaceAlpha',...
        0.4, 'EdgeAlpha', 0.4);
    
    h2 =  histogram(mtlocs,0:5:250,'FaceColor', mtcolor);
    line([median(mtlocs); median(mtlocs)], [0; 25], 'Color', mtcolor, 'LineStyle', '--',...
        'LineWidth', loclinewidth)
    line([median(wtlocs); median(wtlocs)], [0; 25], 'Color', wtcolor, 'LineStyle', ':',...
        'LineWidth', loclinewidth)
    
    hold off
    legend([h1, h2], {'Control', genotype})
end
xlim([15 180])
title('Inter-Spike Interval');
ylabel('# of Events');
xlabel('Interval (s)');

t2.Title.String =['\it' genotype];
t2.TileSpacing = 'Compact';
t2.Padding = 'Compact';

nexttile([1 1])
plotEchos(wormdata,peakthreshold)
title('Inter-Interval Correlation')

nexttile([1 1])
wtcv = plotCV(wtdata, peakthreshold, 0);
mtcv = plotCV(wormdata, peakthreshold, 0);




if iswildtype == 1
    boxplot(wtcv)
    hold on
    scatter(ones(length(wtcv),1),wtcv)
    hold off
    ylim([0 100])
    ax = gca;
    ax.XTickLabel = {'Control'};
    title('Coefficient of Variation')
    ylabel('% of mean')


elseif iswildtype == 0
    cvmat = NaN(max(length(wtcv),length(mtcv)),2);
    cvmat(1:length(wtcv),1) = wtcv;
    cvmat(1:length(mtcv),2) = mtcv; 

    boxplot(cvmat)
    hold on
    scatter(ones(length(cvmat),1),cvmat(:,1))
    scatter(repmat(2,length(cvmat),1),cvmat(:,2))
    hold off


    ylim([0 100])
    ax = gca;
    ax.XTickLabel = {'Control', wormdata(1).genotype};
    title('Coefficient of Variation')
    ylabel('% of mean')
end

end
function shutOffLabels(ax)
set(ax, 'XTickLabels', [], 'YTickLabels', [], 'Box', 'off','Clipping', 'off')
end




function plotGroupAxialMean(axmean, savepath, genotypes,ylimits,axialWindow)
axsavename = savepath;
figure('Position', [89 44.2000 1.2176e+03 688.8000])
tt = tiledlayout('flow');
for a = 1:length(axmean)
    nexttile()
    imagesc(axmean{a}',ylimits)
    title(genotypes{a})
    
    ca = gca;
    
    
ca.XTick = linspace(1,size(axmean{a},1), 7);
ca.XTickLabels = floor(linspace(-axialWindow,axialWindow,7));
    xlabel('Time (s)')
    
    ca.YTick = [ca.YTick(1), ca.YTick(end)];
    ca.YTickLabels = {'Ant','Post'};
    
    
end
tt.Title.String ='Mean Axial Ca^2^+ Signal';
exportgraphics(gcf, axsavename, 'Resolution', 300);
end

function plotGroupCI(CI, CIsavename, genotypes)
figure('Position', [89 44.2000 1.2576e+03 688.8000])
t = tiledlayout('flow');
for i = 1: length(genotypes)
    
    thisCI = CI{i};
    thisGenotype = genotypes{i};

    
    if ~isempty(thisCI)
        
        nexttile();
        scatter(thisCI(:,1), thisCI(:,2))
        %     lsline()
        xlim([0 150])
        ylim([1*10^6 4*10^6])
        title(thisGenotype);
        ylabel('Cumulative Intensity (au)');
        xlabel('Time Since Previous Spike (s)');
    end
end
exportgraphics(gcf, CIsavename, 'Resolution', 300);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plotGroupISI(CI, ISIsavename, genotypes, wtIdx)
figure('Position', [89 44.2000 1.2576e+03 688.8000])
t = tiledlayout('flow');

wtCI = CI{wtIdx};
wtLocs = wtCI(:,1);
for i = 1: length(genotypes)
    thisCI = CI{i};
    thisLocs = thisCI(:,1);
    
    iswildtype = 0;
    thisGenotype = genotypes{i};
    loclinewidth = 1.5;
    wtcolor = [0 0 0];
    mtcolor = [0.09 0.35 0.92];
    if i == wtIdx
        iswildtype =1;
    end
nexttile();
if iswildtype == 1
    histogram(thisLocs,0:5:250,'FaceColor', wtcolor);
    %     ylim([0 25])
    %     xlim([0 250])
    line([median(thisLocs); median(thisLocs)], [0; 25], 'Color', wtcolor, 'LineStyle', ':',...
        'LineWidth', loclinewidth)
    
    
elseif iswildtype == 0
    hold on
    h1=histogram(wtLocs,0:5:250,'FaceColor', [0.7 0.7 0.7],'FaceAlpha',...
        0.4, 'EdgeAlpha', 0.4);
    
    h2 =  histogram(thisLocs,0:5:250,'FaceColor', mtcolor);
    line([median(thisLocs); median(thisLocs)], [0; 25], 'Color', mtcolor, 'LineStyle', '--',...
        'LineWidth', loclinewidth)
    line([median(wtLocs); median(wtLocs)], [0; 25], 'Color', wtcolor, 'LineStyle', ':',...
        'LineWidth', loclinewidth)
    
    hold off
    legend([h1, h2], {'Control', thisGenotype})
end
xlim([15 180])
title(['\it' thisGenotype]);
ylabel('# of Events');
xlabel('Interval (s)');

end
exportgraphics(gcf, ISIsavename, 'Resolution', 300);
end
function plotGroupPR(groupPR, PRsavename, genotypes, axialWindow)
figure('Position', [89 44.2000 1.2576e+03 688.8000])
t = tiledlayout('flow');
for i = 1: length(genotypes)
    
    thisPR = groupPR{i};
    thisGenotype = genotypes{i};
    if ~isempty(thisPR)
        
        nexttile();
        histogram(thisPR,(axialWindow/2)*-1:0.5:axialWindow/2)
        title(thisGenotype);
        ylim([0 20]);
        ylabel('Count');
        xlabel('Length/s (from tail to head)');
    end
end
exportgraphics(gcf, PRsavename, 'Resolution', 300);
end

function plotGroupTraces(groupTraces, GTsavename, genotypes, traceylimits, wtIdx)
singlewidth = 0.5;
tracewidth = 2;
tracelinewidth = 1.5;
loclinewidth = 1.5;
  wtcolor = [0 0 0];
    mtcolor = [0.09 0.35 0.92];
  

wtprofiles = groupTraces{wtIdx};
pktime = linspace(-15,15,length(wtprofiles))';
pktimemat = repmat(pktime, [1,size(wtprofiles,2)]);
wtpk = max(mean(wtprofiles,2,'omitnan'));





figure('Position', [89 44.2000 1.2576e+03 688.8000])
t = tiledlayout('flow');


for i = 1:length(genotypes)
    
mtprofiles = groupTraces{i};
mtpktime = linspace(-15,15,length(mtprofiles))';
mtpktimemat = repmat(mtpktime, [1,size(mtprofiles,2)]);
mtpk = max(mean(mtprofiles,2,'omitnan'));    

if i == wtIdx
    iswildtype =1; 
else 
    iswildtype = 0;
end    
    
nexttile([1 1])


hold on
plot(pktimemat, wtprofiles, 'Color', [0.6 0.6 0.6], 'LineWidth', singlewidth) % plot wild type individual traces
line([-15; 15], [wtpk; wtpk], 'Color', wtcolor, 'LineStyle', ':',...    %line at wild type mean
    'LineWidth', tracelinewidth)

if iswildtype==0
    plot(mtpktimemat, mtprofiles, 'Color', [0.7,0.85,0.99],'LineWidth', singlewidth); % plot mutant individual traces
    line([-15; 15], [mtpk; mtpk], 'Color', mtcolor, 'LineStyle', '--',... % line at mutant mean
        'LineWidth', tracelinewidth)
end

plot(pktime, mean(wtprofiles,2,'omitnan'), 'Color', wtcolor, 'LineWidth', tracewidth); % plot wild type trace mean
if iswildtype ==0
    plot(mtpktime, mean(mtprofiles,2,'omitnan'), 'Color', mtcolor, 'LineWidth', tracewidth); % plot mutant trace mean
end
thisGenotype = genotypes{i};

xlim([-15 15])
ylim(traceylimits)
title(['\it' thisGenotype]);
ylabel('Signal Intensity (a.u.)')
xlabel('Time (s)');
hold off
end
exportgraphics(gcf, GTsavename, 'Resolution', 300);
end

