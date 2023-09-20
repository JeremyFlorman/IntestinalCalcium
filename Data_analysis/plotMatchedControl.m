
parentfolder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_Mutants\slcf-1';


dd = dir(parentfolder);
ddirflag = [dd.isdir];
dd = dd(ddirflag);
dd = dd(3:end);

if nnz(ddirflag) == 2 % see if there are no subfolders
    dd = 1;
end

for q = 1:length(dd)

    if isstruct(dd)     % if there are subfolders, loop through them
    datafolder = [dd(q).folder '\' dd(q).name];
    elseif dd == 1      % if no subfolders, just analyze the parentfolder
        datafolder = parentfolder;
    end
%     datafolder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\egl-19(lf)';
settings = returnPlotSettings();

peakthreshold = settings.peakthreshold; 
traceYLimit = settings.traceylimit; 
axylimit = settings.axylimit;
overlayplots = 0; 


tolimit = 10;             % set to 0 if you want to plot all bulk & axial signal plots.
                         %  set to -1 if you want equal # of control and
                         %  mutant plots, the latter is better for
                         %  comparison as bulk signals will haveS identical
                         %  y-axis scaling. To plot a specific number of
                         %  plots, set tolimit to that number
                      

d = dir([datafolder '\*_mergedData.mat']);
gn = regexp(datafolder, '\', 'split'); % find the mutant genotype name
if ~strcmp(gn{end}, '') 
    mtname = gn{end};  
elseif ~strcmp(gn{end-1}, '')
    mtname = gn{end-1};
end
    
datanames = {d(:).name};
mtflag = contains(datanames, mtname); 
 
if nnz(mtflag) == 1
    mtpath = fullfile(d(mtflag).folder, d(mtflag).name); % find the mutant genotype path
end 


[mtdata, wtdata] = parseWormData(mtpath);


if strcmp(mtdata(1).genotype,'wildtype')== 1
    plotcontrol = 0; % dont plot additional control data if this is the control dataset
    spikeprops = getSpikeLocs(mtdata,[],0);
    mtdata = mtdata(spikeprops.mtsort);

else 
    plotcontrol = 1; % if this is not the control dataset include matched control data. 
end

    if tolimit ==  -1
        if plotcontrol == 1
        plotlimit = min(length(wtdata),length(mtdata));
        elseif plotcontrol == 0
            plotlimit = 4; % how many control plots if there is no mutant data?
        end
    elseif tolimit == 0 
        plotlimit = length(mtdata);
    elseif tolimit > 0
        plotlimit = tolimit;
    end



%% 

 figure('Position', [31.4000 113 1300 662.4000],'Color',[1 1 1])
% figure('Position', [497 88.2000 513.6000 662.4000],'Color',[1 1 1])

t = tiledlayout(5,4, "TileSpacing","compact",  "Padding","tight");

if plotcontrol == 0 % add traces for control matched control data if this is the control dataset
    if plotlimit >0 && plotlimit*2<length(mtdata)
        chunk1 = mtdata(1:plotlimit);
        chunk2 = mtdata(plotlimit+1:plotlimit*2);
    else
        chunk1 = mtdata(1:floor(length(mtdata)/2));
        chunk2 = mtdata(floor(length(mtdata)/2+1):end);
    end
if overlayplots == 0 
nexttile([3 1])
plotBulkSignal(chunk1,0,plotlimit)
 
nexttile([3 1])
plotAxialSignal(chunk1,0,plotlimit)

nexttile([3 1])
plotBulkSignal(chunk2,0,plotlimit)
 
nexttile([3 1])
plotAxialSignal(chunk2,0,plotlimit)

elseif overlayplots == 1
        nexttile([3 2])
        plotOverlay(chunk1,0,plotlimit)
        
        nexttile([3 2])
        plotOverlay(chunk2,0,plotlimit)
end


elseif plotcontrol == 1

    if overlayplots == 0
        nexttile([3 1])
        plotBulkSignal(mtdata, 1,plotlimit)


        nexttile([3 1])
        plotAxialSignal(mtdata, 1,plotlimit)

        nexttile([3 1])
        plotBulkSignal(mtdata, 0,plotlimit)


        nexttile([3 1])
        plotAxialSignal(mtdata, 0,plotlimit)

    elseif overlayplots ==1
        nexttile([3 2])
        plotOverlay(mtdata,1,plotlimit)

        nexttile([3 2])
        plotOverlay(mtdata,0,plotlimit)
    end

end



%% Spike Profiles
nexttile([2 1])
plotSpikeProfiles(mtdata,plotcontrol,1,1)




%% interval histogram 
nexttile([2,1]) 
plotIntervalHistogram(mtdata, plotcontrol,1,1,1)



%% coefficient of variation 
nexttile([2,1])
plotCV(mtdata, plotcontrol,1);

%% interval correlation 
nexttile([2,1])
plotEchos(mtdata, plotcontrol,1,1)

%%
% savename = [datafolder '\' mtname '_Matched_Control_Data.png'];
% savename = strrep(mtpath, 'mergedData.mat', 'Matched_Control_Data.png')
savename = [parentfolder '\' mtname '_matched_Control_Data.png']
exportgraphics(gcf, savename, 'Resolution', 300);
end
% close all