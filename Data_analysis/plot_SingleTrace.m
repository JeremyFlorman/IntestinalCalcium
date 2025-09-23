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
traceylimits= settings.traceylimit;
axylimits= settings.axylimit;
singlespike =settings.singleSpike;
windowInSeconds = settings.spikeWindow;
axSigCMap = settings.axSigCMap;
expLengthInMin = (length(wormdata.bulkSignal)/settings.framerate)/60;
window = settings.framerate*windowInSeconds;
time = linspace(0,expLengthInMin,length(wormdata.bulkSignal));
tracecolor = [0.2 0.2 0.2];
axsig = smoothdata(wormdata.autoAxialSignal,1, 'gaussian', 65);
bulkSignal = wormdata.bulkSignal;
loc = wormdata.peakLoc;
pk = wormdata.peakAmplitude;
locinmin = loc/settings.framerate/60;


% if ~isempty(loc)
%     if loc(singlespike)-window >= 1
%         pre = loc(singlespike)-window;
%         post = loc(singlespike)+window;
%     elseif loc(singlespike)-window < 1
%         singlespike = singlespike+1;
%         pre = loc(singlespike)-window;
%         post = loc(singlespike)+window;
%     end
% else 
    pre = window;
    post = window*2;
% end

for idx = singlespike %:length(loc)
    singlespike = idx;
    figure('Position', [49.8000 222.6000 1060 407.2000], 'Color', [1 1 1])
    t = tiledlayout(2,7,'Padding','tight');
    nexttile([1 5])

    plot(time,bulkSignal, 'Color', tracecolor, 'LineWidth', 1)

    set(gca, 'Box', 'off','Clipping', 'off')
    xlim([0 expLengthInMin])
    xt = get(gca, 'XTick');
    xt = xt*settings.framerate*60;
    xtl = get(gca,'XTickLabels');
    if ~isempty(pk)
        hold on
        plot(locinmin, pk*1.15, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.2 0.2 0.2]);
        plot(locinmin(singlespike), pk(singlespike)*1.15, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.8 0.2 0.2]);
        hold off
    end
    ylim(traceylimits)


    normalization = settings.normalize;
    titles = {'Mean Fluorescence (a.u.)','\Delta F/F (a.u)', 'Z-Score (s.d.)', 'Normalized Fluorescence (a.u.)'};
    if strcmp(normalization, 'None') == 1
        ylabel(titles(1))
    elseif strcmp(normalization, 'Delta F/F0') == 1
        ylabel(titles(2))
    elseif strcmp(normalization, 'Z-Score') == 1
        ylabel(titles(3))
    elseif strcmp(normalization, 'Control') == 1
        ylabel(titles(4))
    end
    xlabel('Time (min)')

    % % % % % Bulk Signal Single spike % % % %


    % % % % % % Axial Signal single spike % % % % %
    nexttile([2,2])
    axsingle = smoothdata(wormdata.autoAxialSignal(pre:post,:),1, 'gaussian', 45);

    colormap(axSigCMap);
    imagesc(axsingle',axylimits)

    set(gca,'XTick',linspace(1,length(axsingle),5),...
        'XTickLabels', linspace(-windowInSeconds,windowInSeconds,5),...
        'YTick', [25 size(axsingle,2)-25], 'YTickLabels', {'Head', 'Tail'});
    xlabel('Time (s)');
    title('axialSignal');
    box off
    cb = colorbar;
    cb.Label.String ='Fluorescent Intensity (a.u.)';
    cb.Label.Rotation = -90;



    % % % % %  % axial signal % % % % % % % linspace(1, length(axsig),length(xtl))
    nexttile([1,5])
    imagesc(axsig', axylimits)
    set(gca, 'XTick',xt,'XTickLabels',xtl,...
        'YTick', [25 size(axsig,2)-25], 'YTickLabels', {'Head', 'Tail'});
    colormap(gca, axSigCMap);
    xlabel('Time (min)')

    if isfield(wormdata, 'pBoc')
        for i =1:length(wormdata.pBoc)
            text(wormdata.pBoc(i), size(axsig, 2)*.05, 'P', 'Color',[1 1 1])
        end
    end

    if isfield(wormdata, 'exp')
        for i =1:length(wormdata.exp)
            text(wormdata.exp(i), size(axsig, 2)*.05, 'x', 'Color',[1 1 1])
        end
    end


    title(t, strrep(name, '_wormdata.mat', ''), 'Interpreter', 'none')

    outpath = [filepath '\' name 'spike_' num2str(singlespike) '_SingleTrace.png'];
    exportgraphics(gcf, outpath)
end

end