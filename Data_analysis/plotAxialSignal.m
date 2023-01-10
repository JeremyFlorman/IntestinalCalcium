function [] = plotAxialSignal(datapath,plotcontrol, plotlimit)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



% datapath = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\flr-1\flr-1_mergedData.mat";
% plotcontrol = 1;
% plotlimit = 3;

[mtdata, wtdata] = parseWormData(datapath);
data = [];

settings = returnPlotSettings();
axylimit = settings.axylimit;
peakthreshold = settings.peakthreshold;

spikeProperties = getSpikeLocs(datapath,peakthreshold,1);
mtsort = spikeProperties.mtsort;
wtsort = spikeProperties.wtsort;

switch plotcontrol
    case 0
        data = mtdata(mtsort);
    case 1
        data = wtdata(wtsort);
end


if plotlimit == 0 || plotlimit>length(data)
    num2plot = length(data);
else
    num2plot = plotlimit;
end


buffer = NaN(length(data(1).autoAxialSignal), 10);

for i = 1:num2plot
    if isfield(data, 'backgroundSignal')
    background = repmat(data(i).backgroundSignal,1,size(data(i).autoAxialSignal,2));
    axsig = data(i).autoAxialSignal-background;
    else
    axsig = data(i).autoAxialSignal;
    end

    if i == 1
        axialMatrix = axsig;
    elseif i>1 
    axialMatrix = horzcat(axialMatrix,buffer, axsig);
    end
end


% figure('Position', [488 70.6000 465 691.4000]);
imagesc(smoothdata(axialMatrix,'movmean',75)',axylimit)
% map = colormap('bone');
% colormap(flipud(map))
% colormap('turbo')
colormap(viridis)

if isfield(data, 'genotype')
    title(['\it' data(1).genotype])
end


ax = gca;
ax.YTickLabel = [];
ax.XTickLabel = num2str(round(str2double(ax.XTickLabel)/15/60,0));
xlabel('Time (min)');
end



