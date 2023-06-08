function [plotSettings] = returnPlotSettings(plotSettings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if nargin == 1
    plotSettings = plotSettings;
else
    plotSettings.peakdistance = 150;
    plotSettings.peakwidth = 15;
    plotSettings.autoFixAxialSignal = 1;
    plotSettings.axSigToQuerry = 0.10;
    plotSettings.framerate = 15;

    plotSettings.axylimit = [-500 12000];
    plotSettings.wtcolor = [0.4 0.4 0.4];
    plotSettings.mtcolor = [0.66 0.74 0.91];
    plotSettings.mtedgecolor = [0.66 0.74 0.91] -.4;
    plotSettings.loclinewidth = 1.5;
    plotSettings.binedges = 0:2:65;     % for histograms/correlation
    plotSettings.xlimits = [0 10];  % for traces
    plotSettings.spikeProfileWindow = 40;
    plotSettings.showFitParams = 1;

    plotSettings.axialXticint = 2; % # x tick interval (in minutes) for axial signal plots.

    plotSettings.sortType = 2; % 0=dont sort,  1=num spikes, 2=amplitude
    plotSettings.sortDir = 'descend';
    plotSettings.normalize = 0;

    plotSettings.tolimit =  10;              % set to 0 if you want to plot all bulk & axial signal plots.
                                             %  set to -1 if you want equal # of control and
                                             %  mutant plots, the latter is better for
                                             %  comparison as bulk signals will haveS identical
                                             %  y-axis scaling. To plot a specific number of
                                             %  plots, set tolimit to that number

     plotSettings.overlayplots = 0; 
     plotSettings.controlname ='wildtype-30mM-5HT';

     plotSettings.numColumns = 1;

    switch plotSettings.normalize
        case 1
            plotSettings.traceylimit = [-.2 1];
            plotSettings.peakthreshold = .1;
        case 0
            plotSettings.traceylimit = [-500 5000];
            plotSettings.peakthreshold = 750;
    end

end
end