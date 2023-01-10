
parentfolder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants';
dd = dir(parentfolder);
ddirflag = [dd.isdir];
dd = dd(ddirflag);
dd = dd(3:end);
% figpos = [77.8000 480.2000 1.2952e+03 247.2000]; % no triple

settings = returnPlotSettings();

peakthreshold = settings.peakthreshold; 
traceYLimit = settings.traceylimit; 
axylimit = settings.axylimit;
plotlimit = 5;
plotnum = [];
plotcontrol = [];


% % all genotypes
% figpos =[74.6000 81.8000 1.3352e+03 660];   % all genotypes
% figpos2 = [20.2000 89 1.3512e+03 672];
% 
% genotypes = {'wildtype', 'cca-1', 'egl-19(gf)','egl-19(lf)', 'trpa-1','flr-1', ...
%              'inx-16', 'eat-2', 'tph-1', 'tdc-1','lgc-55;ser-2;tyra-2;tyra-3','egl-8', ...
%             'plc-3(kd)','itr-1','unc-43(n498sd)', 'unc-43(sa200)', 'unc-43(e408)', 'itr-1-unc-43',...
%             'dec-1','dec-2', 'dec-7','dec-9','dec-10',};
% w = 6;
% h = ceil(length(genotypes)/w);
% prefix = '\all_Genotypes ';

% fewer genotypes
figpos = [79.4000 568.2000 1.3784e+03 210.4000];
figpos2 = [130.6000 345.8000 1324 325.6000];

%  genotypes = {'wildtype','unc-43(e408)','itr-1-unc-43', 'itr-1'};
% genotypes = {'wildtype', 'dec-7', 'dec-9','dec-10'};
genotypes = {'wildtype', 'dec-9', 'unc-43(sa200)'};

if length(genotypes) == 4
    figpos = [245.8000 117 948 210.4000];
    figpos2 = [305.8000 308.2000 848 325.6000];
elseif length(genotypes) == 3
    figpos = [263.4000 369 761.6000 231.2000];
    figpos2 = [447.4000 320.2000 703.2000 325.6000];
end

w = length(genotypes);
h = ceil(length(genotypes)/w);
prefix = '\dec_mutants ';

%%

plotorder = nan(length(genotypes),1);
for i = 1:length(genotypes)
  plotorder(i) = find(strcmp({dd.name},genotypes{i}));
end




spikefig = figure('Position', figpos);
spiket = tiledlayout(h, w,'Parent', spikefig, 'TileSpacing','compact','Padding','compact');

histfig = figure('Position', figpos);
histt = tiledlayout(h, w,'Parent', histfig, 'TileSpacing','compact','Padding','compact');

cvfig = figure('Position', figpos);
cvt = tiledlayout(h, w,'Parent', cvfig, 'TileSpacing','compact','Padding','compact');

correlationfig = figure('Position', figpos);
cort = tiledlayout(h, w,'Parent', correlationfig, 'TileSpacing','compact','Padding','compact');

bulkfig = figure('Position', figpos2, 'Color', [1 1 1]);
bulkt = tiledlayout(h, w,'Parent', bulkfig, 'TileSpacing','tight','Padding','tight');

axialfig = figure('Position', figpos2, 'Color', [1 1 1]);
axt = tiledlayout(h, w,'Parent', axialfig, 'TileSpacing','tight','Padding','tight');

overlayfig = figure('Position', figpos2, 'Color', [1 1 1]);
olt = tiledlayout(h, w,'Parent', overlayfig, 'TileSpacing','tight','Padding','tight');




%  plotorder = [7 2 6 4 1]; % tyramine no triple mutants
% plotorder = [7 2 6 1 3 5 4]; % tyramine
% plotorder = [16 1 4 5 12 6 3 11 8 7 9 14 15 13 10 2]; %dmp mutants.





for q = 1:length(plotorder)
   

    plotnum = plotorder(q);   % this specifies what order to plot genotypes... 
    datafolder = [dd(plotnum).folder '\' dd(plotnum).name]; %change plotnum to q to plot in default order. 


        if  mod(q-1,w) == 0
            labelYAxis = 1;
        else 
            labelYAxis = 0;
        end
        
        if q > h*w-w
            labelXAxis = 1;
        else 
            labelXAxis = 0;
        end

d = dir([datafolder '\*_mergedData.mat']);
gn = regexp(datafolder, '\', 'split'); % find the mutant genotype name
mtname = gn{end};   % find the mutant genotype name
datanames = {d(:).name};
mtflag = contains(datanames, mtname); 

if nnz(mtflag) == 1
    mtpath = fullfile(d(mtflag).folder, d(mtflag).name); % find the mutant genotype path
end
 
[mtdata, wtdata] = parseWormData(mtpath);

% 
if strcmp(mtdata(1).genotype, 'wildtype') ==1 % make sure the first plot is wildtype
    plotcontrol = 0;
else
    plotcontrol = 1;
end


% if plotnum == 1
%     plotcontrol = 1;
% else 
%     plotcontrol = 0;
% end

%% spike profiles
nexttile(spiket, q)
plotSpikeProfiles(mtdata,plotcontrol, labelXAxis,labelYAxis)
title(['\it' mtdata(1).genotype])




%% interval histogram 
nexttile(histt, q)  
plotIntervalHistogram(mtdata, plotcontrol,0,labelXAxis,labelYAxis)
title(['\it' mtdata(1).genotype])


%% coefficient of variation 
nexttile(cvt,q)
plotCV(mtdata, plotcontrol,labelYAxis);
title(['\it' mtdata(1).genotype])

%% interval correlation 
nexttile(cort,q)
plotEchos(mtdata, 0,labelXAxis,labelYAxis)
title(['\it' mtdata(1).genotype])



    %% axial signal 
    
 
    nexttile(axt, q) 
    plotAxialSignal(mtdata,0,plotlimit)
    title(['\it' mtdata(1).genotype],'FontSize', 8)
    
    ax2 = gca;
    ax2.YTick = [];
    ax2.XTick = [];
    ax2.XLabel = [];

%% bulk signal

    nexttile(bulkt, q)
    plotBulkSignal(mtdata,0,plotlimit)
     title(['\it' mtdata(1).genotype],'FontSize', 8)
    ax = gca;
    ax.YTick = []; 
    ax.XTick = [];
    ax.XLabel = [];
%% Overlay Signal
    nexttile(olt, q)
    plotOverlay(mtdata,0,plotlimit)
    title(['\it' mtdata(1).genotype],'FontSize', 8)
    ax0 = gca;
    ax0.YTick = []; 
    ax0.XTick = [];
    ax0.XLabel = [];

end


% 
% title(spiket, 'Spike Profiles')
% title(histt,'Inter-spike Interval');
% title(cvt, 'Coefficient of Variation');
% title(cort,'Inter-interval Correlation');
% title(bulkt,'Bulk Ca^2^+ Signal');
% title(axt,'Axial Ca^2^+ Signal');


exportgraphics(spikefig, [parentfolder prefix 'Spike_Profiles.png'], 'Resolution', 300);
exportgraphics(histfig, [parentfolder prefix 'Interval_Histogram.png'], 'Resolution', 300);
exportgraphics(cvfig, [parentfolder prefix 'Coefficent_of_Variance.png'], 'Resolution', 300);
exportgraphics(correlationfig, [parentfolder prefix 'Interval_Correlation.png'], 'Resolution', 300);
exportgraphics(bulkfig, [parentfolder prefix 'Bulk_Signal.png'], 'Resolution', 300);
exportgraphics(axialfig, [parentfolder prefix 'axial_Signal.png'], 'Resolution', 300);
exportgraphics(overlayfig, [parentfolder prefix 'Overlay_Signal.png'], 'Resolution', 300);

