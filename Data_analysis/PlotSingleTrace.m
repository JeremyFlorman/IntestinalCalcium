[file ,path] = uigetfile('Y:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\Receptor_Mutants\wildtype-30mM-TA\*.mat')

% file = '220316_zfis178_itr-1_4_wormdata.mat';
% path = 'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\itr-1\220316_zfis178_itr-1_4';
filename = [path '\' file];
[fp, name, ext] = fileparts(filename);

singlespike =5;
window = 60*15;

tempdir = 'C:\tmp';


wormdata = copyLoadClear(filename, tempdir);
wormdata = wormdata.wormdata;
%%
settings = returnPlotSettings();

peakthreshold=settings.peakthreshold;
peakdistance = settings.peakdistance;
traceylimits= settings.traceylimit;
axylimits= settings.axylimit;
% axylimits=[0 35000];

time = linspace(0,10,length(wormdata.bulkSignal));

tracecolor = [0.2 0.2 0.2];


if isfield(wormdata, 'backgroundSignal')
    background = wormdata.backgroundSignal;
    rawsig = wormdata.bulkSignal;

    sig = rawsig-background;
else
    sig = wormdata.bulkSignal;
end




terpdata = fillmissing(sig, 'movmedian',100);
[pk, loc] = findpeaks(terpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance);
locinmin = loc/15/60;

for idx = singlespike %:length(loc)
    singlespike = idx;
    figure('Position', [81.8000 477.8000 1.0776e+03 284.2000], 'Color', [1 1 1])
    t = tiledlayout(2,5,'Padding','tight');
    nexttile([1 3])

    plot(time,terpdata, 'Color', tracecolor, 'LineWidth', 1)
    xtl = get(gca,'XTickLabels');
    set(gca, 'XTickLabels', [], 'Box', 'off','Clipping', 'off')

    xlim([0 10])
    if ~isempty(pk)
        hold on
        plot(locinmin, pk*1.15, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.2 0.2 0.2]);
        plot(locinmin(singlespike), pk(singlespike)*1.15, 'v', 'MarkerSize', 6, 'MarkerEdgeColor', [0.8 0.2 0.2]);
        hold off
    end
    %     ylim(traceylimits)
    ylabel('Signal Intensity (a.u.)')
%     title('Bulk Signal')

    %%%%%%%%% axial signal
    if ~isfield(wormdata, 'autoAxialSignal')
        if isfield(wormdata, 'rawAxialSignal')
            afs = autoFixSignal(wormdata.rawAxialSignal);
            wormdata.autoAxialSignal = afs;
        end
    end

    axsig = wormdata.autoAxialSignal;
    if isfield(wormdata, 'backgroundsignal')
    background = repmat(wormdata.backgroundSignal,1,size(axsig,2));
    axsig = axsig-background;
    end

    %% single spike
    nexttile([2,2])
    
    if ~isempty(loc)
    pre = loc(singlespike)-1100%window;
    post = loc(singlespike)+750%window;
    else
        pre = 900-window;
        post = 900+window;
    end
    axsingle = axsig(pre:post,:);
    %         axsingle(49:52,:) = fliplr(axsingle(49:52,:));
    imagesc(smoothdata(axsingle,'movmean',15)',axylimits)
    xt = get(gca, 'XTick');
    xlab  = {};
    for idx = 1:length(xt)
        n = num2str(floor(xt(idx)/15));
        xlab(idx) = {n};
    end
    %         xlab = floor(linspace(1,(xt(end)-xt(1))/15,length(xt)));
    set(gca, 'XTicklabels',xlab, 'YTick', [25 size(axsingle,2)-25], 'YTickLabels', {'Head', 'Tail'},'FontSize',9, 'FontWeight','bold');
    xlabel('Time (s)');
    title('Single Spike');
    colormap(gca, 'turbo')

    
    
    %% axial signal
    nexttile([1,3])
    imagesc(smoothdata(axsig,'movmean',75)', axylimits)
    set(gca, 'XTickLabels', round(get(gca,'XTick')/900), 'YTick', [25 size(axsig,2)-25], 'YTickLabels', {'Head', 'Tail'},'FontSize',9, 'FontWeight','bold');
    colormap(gca, 'bone');
    xlabel('Time (min)')

%     title('Axial Signal');



    title(t, strrep(file, '_wormdata.mat', ''), 'Interpreter', 'none')

    outpath = [path '\' name 'spike_' num2str(singlespike) '_SingleTrace.png'];
    exportgraphics(gcf, outpath)
end