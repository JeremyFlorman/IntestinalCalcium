%  [file ,path] = uigetfile('Y:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\Receptor_Mutants\wildtype-30mM-TA\*.mat')

% filename = 'Y:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\Receptor_Mutants\wildtype-30mM-TA\211108_zfis178_wildtype-control_2\211108_zfis178_wildtype-control_2_wormdata.mat';
% filename = 'Y:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\Receptor_Mutants\wildtype-30mM-TA\220126_zfis178_wildtype-30mM-TA_1\220126_zfis178_wildtype-30mM-TA_1_wormdata.mat';
filename = 'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\slcf-1\230505_zfis178_wildtype_4\230505_zfis178_wildtype_4_wormdata.mat'
[path, name, ext] = fileparts(filename);
file = [name ext];



tempdir = 'C:\tmp';


wormdata = copyLoadClear(fullfile(path,file), tempdir);
wormdata = wormdata.wormdata;
%%
settings = returnPlotSettings();

peakthreshold=settings.peakthreshold;
peakdistance = settings.peakdistance;
traceylimits= settings.traceylimit;
axylimits= settings.axylimit;
traceylimits = [0 8000];
% axylimits=[0 20000];


singlespike =3;
windowInSeconds = 10;
window = 15*windowInSeconds;

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
    figure('Position', [49.8000 222.6000 1060 407.2000], 'Color', [1 1 1])
    t = tiledlayout(2,7,'Padding','tight');
    nexttile([1 3])

    plot(time,terpdata, 'Color', tracecolor, 'LineWidth', 1)
    xtl = get(gca,'XTickLabels');
    
    set(gca, 'Box', 'off','Clipping', 'off')

    xlim([0 10])
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

    %%%%%%%%% axial signal
    if ~isfield(wormdata, 'autoAxialSignal')
        if isfield(wormdata, 'rawAxialSignal')
            afs = autoFixSignal(wormdata.rawAxialSignal);
            wormdata.autoAxialSignal = afs;
        end
    end 

    axsig = smoothdata(wormdata.autoAxialSignal, 'movmean',30);
    if isfield(wormdata, 'backgroundSignal')
    backgroundMatrix = repmat(wormdata.backgroundSignal,1,size(axsig,2));
    
    axsig = axsig-backgroundMatrix;
    end

    % % % % % Bulk Signal Single spike % % % % 
        
    if ~isempty(loc)
    pre = loc(singlespike)-window;
    post = loc(singlespike)+window;
    else
        pre = 900-window;
        post = 900+window;
    end
    
    nexttile([2,2])
    singleBulk = terpdata(pre:post,1);
    [spk, sloc] = findpeaks(singleBulk, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance);

    singleTime = linspace(-windowInSeconds,windowInSeconds,length(singleBulk));
    plot(singleTime,singleBulk, 'Color', tracecolor, 'LineWidth', 1)
    hold on
    if ~isempty(pk)
    plot(0, pk(singlespike)*1.05, 'v', 'MarkerSize', 10, 'MarkerEdgeColor', [0.2 0.2 0.2])
    hold off
    end
    
    ax = gca;
    ax.XTick = linspace(-windowInSeconds, windowInSeconds, 5);
    ax.XTickLabel = linspace(-windowInSeconds, windowInSeconds, 5);
    ax.XLim = [singleTime(1) singleTime(end)];
    title('Bulk Signal')
    xlabel('Time (s)')
    ylabel('Signal Intensity (a.u.)');
    box off
    
    % % % % % % Axial Signal single spike % % % % % 
    nexttile([2,2])
    axsingle = axsig(pre:post,:);
    
    colormap('viridis')
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
    nexttile([1,3])
    imagesc(axsig', axylimits)
    set(gca, 'XTick', linspace(1, length(axsig),length(xtl)),'XTickLabels',xtl,...
        'YTick', [25 size(axsig,2)-25], 'YTickLabels', {'Head', 'Tail'});
    colormap(gca, 'viridis');
    xlabel('Time (min)')

%     title('Axial Signal');



    title(t, strrep(file, '_wormdata.mat', ''), 'Interpreter', 'none')

    outpath = [path '\' file 'spike_' num2str(singlespike) '_SingleTrace.png'];
    exportgraphics(gcf, outpath)
end