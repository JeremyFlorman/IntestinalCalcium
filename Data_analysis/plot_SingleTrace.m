function [outputArg1,outputArg2] = plot_SingleTrace(filename,settings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<1
filename = 'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\slcf-1\230505_zfis178_wildtype_4\230505_zfis178_wildtype_4_wormdata.mat';
end

if nargin<2
    settings = returnPlotSettings();
end

[filepath,name,~] = fileparts(filename);

tempdir = 'C:\tmp';


wormdata = copyLoadClear(filename, tempdir);
wormdata = wormdata.wormdata;

if length(wormdata) >1
    wormdata = wormdata(4);
end

[wormdata, ~, ~] = processWormdata(wormdata,settings);
%%
% settings = returnPlotSettings();

% peakthreshold=settings.peakthreshold;
% peakdistance = settings.peakdistance;
traceylimits= settings.traceylimit;
axylimits= settings.axylimit;
singlespike =settings.singleSpike;
windowInSeconds = settings.spikeWindow;
axSigCMap = settings.axSigCMap;
expLengthInMin = ceil((length(wormdata.bulkSignal)/settings.framerate)/60);
% traceylimits = [0 8000];
% axylimits=[0 20000];
 


window = settings.framerate*windowInSeconds;

time = linspace(0,expLengthInMin,length(wormdata.bulkSignal));

tracecolor = [0.2 0.2 0.2];


% if isfield(wormdata, 'backgroundSignal')
%     background = wormdata.backgroundSignal;
%     rawsig = wormdata.bulkSignal;
% 
%     sig = rawsig-background;
% else
%     sig = wormdata.bulkSignal;
% end
% 
% 

axsig = smoothdata(wormdata.autoAxialSignal,1, 'gaussian', 30);
bulkSignal = wormdata.bulkSignal;
loc = wormdata.peakLoc;
pk = wormdata.peakAmplitude;
locinmin = loc/15/60;

for idx = singlespike %:length(loc)
    singlespike = idx;
    figure('Position', [49.8000 222.6000 1060 407.2000], 'Color', [1 1 1])
    t = tiledlayout(2,7,'Padding','tight');
    nexttile([1 5])

    plot(time,bulkSignal, 'Color', tracecolor, 'LineWidth', 1)
    xtl = get(gca,'XTickLabels');
    
    set(gca, 'Box', 'off','Clipping', 'off')

    xlim([0 expLengthInMin])
    if ~isempty(pk)
        hold on
        plot(locinmin, pk*1.15, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.2 0.2 0.2]);
        plot(locinmin(singlespike), pk(singlespike)*1.15, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.8 0.2 0.2]);
        hold off
    end
    ylim(traceylimits)
    ylabel('Signal Intensity (a.u.)')
    xlabel('Time (min)')
%     title('Bulk Signal')

    % % % % % Bulk Signal Single spike % % % % 
        
    if ~isempty(loc)
    pre = loc(singlespike)-window;
    post = loc(singlespike)+window;
    else
        pre = 900-window;
        post = 900+window;
    end
    
%     nexttile([2,2])
%     singleBulk = bulkSignal(pre:post,1);
% %     [spk, sloc] = findpeaks(singleBulk, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance);
% 
%     singleTime = linspace(-windowInSeconds,windowInSeconds,length(singleBulk));
%     plot(singleTime,singleBulk, 'Color', tracecolor, 'LineWidth', 1)
%     hold on
%     if ~isempty(pk)
%     plot(0, pk(singlespike)*1.05, 'v', 'MarkerSize', 10, 'MarkerEdgeColor', [0.2 0.2 0.2])
%     hold off
%     end
% 
%     ax = gca;
%     ax.XTick = linspace(-windowInSeconds, windowInSeconds, 5);
%     ax.XTickLabel = linspace(-windowInSeconds, windowInSeconds, 5);
%     ax.XLim = [singleTime(1) singleTime(end)];
%     title('Bulk Signal')
%     xlabel('Time (s)')
%     ylabel('Signal Intensity (a.u.)');
%     box off
    
    % % % % % % Axial Signal single spike % % % % % 
    nexttile([2,2])
    axsingle = smoothdata(wormdata.autoAxialSignal(pre:post,:),1, 'gaussian', 15);
    
    colormap(axSigCMap);
    imagesc(axsingle',axylimits)
    xt = get(gca, 'XTick');
    xlab  = {};


    set(gca,'XTick',linspace(1,length(axsingle),5),...
        'XTickLabels', linspace(-windowInSeconds,windowInSeconds,5),...
        'YTick', [25 size(axsingle,2)-25], 'YTickLabels', {'Head', 'Tail'});
    xlabel('Time (s)');
    title('axialSignal');
    box off
    cb = colorbar;
    cb.Label.String ='Fluorescent Intensity (a.u.)'; 
    cb.Label.Rotation = -90;    
       


% % % % %  % axial signal % % % % % % % 
    nexttile([1,5])
    imagesc(axsig', axylimits)
    set(gca, 'XTick', linspace(1, length(axsig),length(xtl)),'XTickLabels',xtl,...
        'YTick', [25 size(axsig,2)-25], 'YTickLabels', {'Head', 'Tail'});
    colormap(gca, axSigCMap);
    xlabel('Time (min)')

%     title('Axial Signal'); 



    title(t, strrep(name, '_wormdata.mat', ''), 'Interpreter', 'none')

    outpath = [filepath '\' name 'spike_' num2str(singlespike) '_SingleTrace.png'];
    exportgraphics(gcf, outpath)
end
end