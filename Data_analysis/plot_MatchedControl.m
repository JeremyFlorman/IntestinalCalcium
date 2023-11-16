function [] = plot_MatchedControl(parentfolder,settings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
parentfolder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\OAS\5-HT\wildtype+5HT';
end

if nargin<2
    settings = returnPlotSettings();
end

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





                     

d = dir([datafolder '\*_mergedData.mat']);
gn = regexp(datafolder, '\', 'split'); % find the mutant genotype name
if ~strcmpi(gn{end}, '') 
    mtname = gn{end};  
elseif ~strcmpi(gn{end-1}, '')
    mtname = gn{end-1};
end
    
datanames = {d(:).name};
mtflag = contains(datanames, mtname); 
 
if nnz(mtflag) == 1
    mtpath = fullfile(d(mtflag).folder, d(mtflag).name); % find the mutant genotype path
elseif length(d)  == 1
    mtpath = fullfile(d.folder, d.name);
end


% [mtdata, wtdata,settings] = parseWormData(mtpath,settings);
[mtdata, wtdata, settings] = processWormdata(mtpath, settings);

tolimit = settings.tolimit;
overlayplots = settings.overlayplots;
controlname = settings.controlname;


if strcmp(mtdata(1).genotype,controlname)== 1
    plotcontrol = 0; % dont plot additional control data if this is the control dataset
    wtdata = [];
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

 figure('Position', [203.4000 86.6000 1.1456e+03 662.4000],'Color',[1 1 1])
% figure('Position', [497 88.2000 513.6000 662.4000],'Color',[1 1 1])

t = tiledlayout(6,4, "TileSpacing","compact","Padding","tight");

if plotcontrol == 0 % add traces for control matched control data if this is the control dataset
    if plotlimit >0 && plotlimit*2<length(mtdata)
        chunk1 = mtdata(1:plotlimit);
        chunk2 = mtdata(plotlimit+1:plotlimit*2);
    else
        chunk1 = mtdata(1:floor(length(mtdata)/2));
        chunk2 = mtdata(floor(length(mtdata)/2+1):end);
    end
if overlayplots == 0 
nexttile([4 1])
plotBulkSignal(chunk1,settings)
 
nexttile([4 1])
plotAxialSignal(chunk1,settings)

nexttile([4 1])
plotBulkSignal(chunk2,settings)
 
nexttile([4 1])
plotAxialSignal(chunk2,settings)

elseif overlayplots == 1
        nexttile([4 2])
        plotOverlay(chunk1,settings)
        
        nexttile([4 2])
        plotOverlay(chunk2,settings)
end


elseif plotcontrol == 1

    if overlayplots == 0
        nexttile([4 1])
        plotBulkSignal(wtdata, settings)


        nexttile([4 1])
        plotAxialSignal(wtdata, settings)

        nexttile([4 1])
        plotBulkSignal(mtdata, settings)


        nexttile([4 1])
        plotAxialSignal(mtdata, settings)

    elseif overlayplots ==1
        nexttile([4 2])
        plotOverlay(wtdata,settings)

        nexttile([4 2])
        plotOverlay(mtdata,settings)
    end

end



%% Spike Profiles
nexttile([2 1])
plotSpikeProfiles(mtdata,wtdata,settings,1,1)




%% interval histogram 
nexttile([2,1]) 
plotIntervalHistogram(mtdata, wtdata, settings,1,1,1)



%% coefficient of variation 
nexttile([2,1])
plotCV(mtdata, wtdata,settings,1);

%% interval correlation  
nexttile([2,1])
plotEchos(mtdata,wtdata,settings,1,1)

%%
% savename = [datafolder '\' mtname '_Matched_Control_Data.png'];
% savename = strrep(mtpath, 'mergedData.mat', 'Matched_Control_Data.png')
savename = [parentfolder '\' mtname '_matched_Control_Data.png']
exportgraphics(gcf, savename, 'Resolution', 300);
end
% assignin('base', 'mtdata', mtdata)
% assignin('base', 'wtdata', wtdata)


end