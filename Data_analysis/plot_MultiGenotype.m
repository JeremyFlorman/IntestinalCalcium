function plot_MultiGenotype(parentfolder, settings)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% prefix = '\dec-mutants'
prefix = settings.prefix;
genotypes = settings.genotypes;
controlname = settings.controlname;
pltOvly = settings.plotOverlay;
pltBulk = settings.plotBulk;
pltAx = settings.plotAxial;
pltCV = settings.plotCV;
pltHist = settings.plotHist;
pltProfile = settings.plotProfile;
pltCorr = settings.plotCorr;

if isempty(genotypes)
mergedDataFiles = dir([parentfolder '\**\*mergedData.mat']);
genotypes = strrep({mergedDataFiles.name},'_mergedData.mat','')';
end


if isfield(settings, 'graphPos')
    figpos = settings.graphPos;
    figpos2 = settings.tracePos;
elseif length(genotypes) == 3
    figpos = [263.4000 369 761.6000 231.2000];
    figpos2 = [771.4000 109.8000 682.4000 585.6000];
elseif length(genotypes) == 4
    figpos = [253.8000 581.8000 600.8000 156];
    figpos2 = [172.2000 233.8000 1.1264e+03 384.8000];

elseif length(genotypes) == 5
    figpos = [135.4000 266.6000 1.0528e+03 241.6000];
    figpos2 = [161 489.8000 1036 174.4000];
else
    figpos = [105 207.4000 1.1672e+03 290.4000];
    figpos2 = [118.6000 181 1.0152e+03 477.6000];
end

w = settings.numColumns;
if w == 0 
    w = length(genotypes);
end

h = ceil(length(genotypes)/w);


%%




if pltProfile
    spikefig = figure('Position', figpos);
    spiket = tiledlayout(h, w,'Parent', spikefig, 'TileSpacing','compact','Padding','compact');
end

if pltHist
    histfig = figure('Position', figpos);
    histt = tiledlayout(h, w,'Parent', histfig, 'TileSpacing','compact','Padding','compact');
end

if pltCV
    cvfig = figure('Position', figpos);
    cvt = tiledlayout(h, w,'Parent', cvfig, 'TileSpacing','compact','Padding','compact');
end

if pltCorr
    correlationfig = figure('Position', figpos);
    cort = tiledlayout(h, w,'Parent', correlationfig, 'TileSpacing','compact','Padding','compact');
end

if pltBulk
    bulkfig = figure('Position', figpos2, 'Color', [1 1 1]);
    bulkt = tiledlayout(h, w,'Parent', bulkfig, 'TileSpacing','compact','Padding','tight');
end

if pltAx
    axialfig = figure('Position', figpos2, 'Color', [1 1 1]);
    axt = tiledlayout(h, w,'Parent', axialfig, 'TileSpacing','compact','Padding','tight');
end

if pltOvly
    overlayfig = figure('Position', figpos2, 'Color', [1 1 1]);
    olt = tiledlayout(h, w,'Parent', overlayfig, 'TileSpacing','compact','Padding','tight');
end



%  plotorder = [7 2 6 4 1]; % tyramine no triple mutants
% plotorder = [7 2 6 1 3 5 4]; % tyramine
% plotorder = [16 1 4 5 12 6 3 11 8 7 9 14 15 13 10 2]; %dmp mutants.





for q = 1:length(genotypes) 

    d = dir([parentfolder '\**\*' genotypes{q} '*_mergedData.mat']);

    if length(d)>1 % if there is more than one file (like in a double mutant) exclude partial matches
        matchIdx = nan(length(d),1);
        for n = 1:length(d)
            matchIdx(n) = strcmpi(strrep(d(n).name,'_mergedData.mat',''), genotypes{q});
        end
        d = d(logical(matchIdx));
    end


    if length(d)>1    % if there is more than one file (like with control) use the biggest.
        [~,idx] = sort([d(:).bytes],'descend');
        d = d(idx(1));
    end
    
    wormDataPath = fullfile(d(1).folder,d(1).name)

    % try
    [mtdata, wtdata, settings] = processWormdata(wormDataPath, settings);
    % catch
        % warning(['Genotype: ' genotypes{q} ' could not be found! Check spelling and OutputDir path.'])
    % end
%%
% datafolder = 'C:\Users\Jeremy\Dropbox\ins-3 intestinal calcium fig\data\forExport';
%     switch q
%         case 1
%             wormdata = {mtdata(1:5).autoAxialSignal};
%             save([datafolder '\wildtypeControlData.mat'], 'wormdata')
%         case 2
%             wormdata = {mtdata(1:5).autoAxialSignal};
%             save([datafolder '\tdc-1ControlData.mat'], 'wormdata');
%         case 3
%             wormdata = {mtdata(1:5).autoAxialSignal};
%             save([datafolder '\tyra-3ControlData.mat'], 'wormdata');
%         % case 4
%         %     wormdata = {mtdata(1:20).autoAxialSignal};
%         %     save([datafolder '\quadData.mat'], 'wormdata');
%     end

  %%  
    %% Write data to text file---> Add to GUI
    writeData = 0;
    if writeData == 1
    % File names
    sortedFileNames = stripPath(mtdata);
    fn =  strrep(fullfile(d(1).folder,d(1).name), '_mergedData.mat', '_sortOrder.txt');
    writecell(sortedFileNames, fn)
    
    % intervals 
    fn =  strrep(fullfile(d(1).folder,d(1).name), '_mergedData.mat', '_intervals.txt');
    wtInt = nan(length(wtdata),1);
    mtInt = nan(length(mtdata),1);
    for i = 1:length(wtdata)
        wtInt(i) = mean(diff(wtdata(i).peakLoc)/settings.framerate,'omitnan');
    end

    for i = 1:length(mtdata)
        mtInt(i) = mean(diff(mtdata(i).peakLoc)/settings.framerate,'omitnan');
    end
    intervals = {wtInt; mtInt};
    writecell(intervals, fn)
    
    % peak traces
    fn =  strrep(fullfile(d(1).folder,d(1).name), '_mergedData.mat', '_peakTraces.xlsx');
    if exist(fn,"file")
        delete(fn)
    end
    writematrix([mtdata(:).peakTraces], fn, sheet=d(1).name)
    writematrix([wtdata(:).peakTraces], fn, sheet="Control");
    end
    %%


    if ~isfield(mtdata,'genotype')
        mtdata(1).genotype = genotypes{q};
    end

    if ~isfield(wtdata, 'genotype')
        wtdata(1).genotype = genotypes{q};
    end

    if strcmp(genotypes(q), controlname) ==1
        wtdata = [];
    end



    %% shut off interior axis labels

    if  mod(q-1,w) == 0
        labelYAxis = 1;
    else
        labelYAxis = 0;
    end

    %     % for multi column figures.
    if q > h*w-w
        labelXAxis = 1;
    else
        labelXAxis = 0;
    end

    % for 1 colum figures

    %     if q<length(genotypes)
    %         labelXAxis = 0;
    %     else
    %         labelXAxis =1;
    %     end

    %% spike profiles
    if pltProfile
        nexttile(spiket, q)
        plotSpikeProfiles(mtdata,wtdata,settings, labelXAxis,labelYAxis)
        title(['\it' mtdata(1).genotype])
    end



    %% interval histogram
    if pltHist
        nexttile(histt, q)
        plotIntervalHistogram(mtdata,wtdata,settings,0,labelXAxis,labelYAxis)
        title(['\it' mtdata(1).genotype])
    end

    %% coefficient of variation
    if pltCV
        nexttile(cvt,q)
        plotCV(mtdata,wtdata,settings,labelYAxis);
        title(['\it' mtdata(1).genotype])
    end

    %% interval correlation
    if pltCorr
        nexttile(cort,q)
        plotEchos(mtdata,[],settings,labelXAxis,labelYAxis)
        
        title(['\it' mtdata(1).genotype])
        daspect([1 1 1])
    end


    %% axial signal
    if pltAx
        nexttile(axt, q)
        plotAxialSignal(mtdata,settings,labelXAxis)
        ttl = strrep(mtdata(1).genotype, '-GFPBeads-noFood','');
        ttl = strrep(ttl, '-GPFBeads-noFood','');
        ttl = strrep(ttl, '-0-5um','');
        title(['\it' ttl],'FontSize', 8)
    end

    %% bulk signal
    if pltBulk
        nexttile(bulkt, q)
        plotBulkSignal(mtdata,settings,labelXAxis)
        ttl = strrep(mtdata(1).genotype, '-GFPBeads-noFood','');
        ttl = strrep(ttl, '-GPFBeads-noFood','');
        ttl = strrep(ttl, '-0-5um','');
        title(['\it' ttl],'FontSize', 8)
    end

    %% Overlay Signal
    if pltOvly
        nexttile(olt, q)
        plotOverlay(mtdata,settings,labelXAxis)
        title(['\it' mtdata(1).genotype],'FontSize', 8)
    end

end


%
% title(spiket, 'Spike Profiles')
% title(histt,'Inter-spike Interval');
% title(cvt, 'Coefficient of Variation');
% title(cort,'Inter-interval Correlation');
% title(bulkt,'Bulk Ca^2^+ Signal');
% title(axt,'Axial Ca^2^+ Signal');

if pltProfile
    exportgraphics(spikefig, [parentfolder prefix 'Spike_Profiles.png'], 'Resolution', 300);
end

if pltHist
    exportgraphics(histfig, [parentfolder '_' prefix '_Interval_Histogram.png'], 'Resolution', 300);
end

if pltCV
    exportgraphics(cvfig, [parentfolder '_' prefix '_Coefficent_of_Variance.png'], 'Resolution', 300);
end

if pltCorr
    exportgraphics(correlationfig, [parentfolder '_' prefix '_Interval_Correlation.png'], 'Resolution', 300);
end

if pltBulk
    exportgraphics(bulkfig, [parentfolder '_' prefix '_Bulk_Signal.png'], 'Resolution', 300);
end

if pltAx
    exportgraphics(axialfig, [parentfolder '_' prefix '_axial_Signal.png'], 'Resolution', 300);
end

if pltOvly
    exportgraphics(overlayfig, [parentfolder '_' prefix '_Overlay_Signal.png'], 'Resolution', 300);
end

end