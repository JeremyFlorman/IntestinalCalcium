function [plotSettings] = returnPlotSettings(plotSettings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if nargin == 1
    plotSettings = plotSettings;
else

    plotsettings.OAS = 1; % are we using openAutoScope? this is an 8bit camera so pixel values need to be adjusted

    plotSettings.peakdistance = 15;
    plotSettings.peakwidth = 15;
    plotSettings.autoFixAxialSignal = 1;
    plotSettings.axSigToQuerry = 0.10;
    plotSettings.framerate = 15;

    plotSettings.axylimit = [0 35000];
    plotSettings.axSigCMap = 'viridis';
    plotSettings.wtcolor = [0.4 0.4 0.4];
    plotSettings.mtcolor = [0.66 0.74 0.91];
    plotSettings.mtedgecolor = [0.66 0.74 0.91] -.4;


    plotSettings.loclinewidth = 1.5;
    plotSettings.binedges = 0:2:90;     % for histograms/correlation
    plotSettings.xlimits = [0 10];  % for traces
    plotSettings.spikeProfileWindow = 40;
    plotSettings.showFitParams = 0;
    plotSettings.validatePropagationRate = 0;
    plotSettings.validateRiseFall = 0;

    plotSettings.axialXticint = 2; % # x tick interval (in minutes) for axial signal plots.

    plotSettings.sortType = 0; % 0=dont sort,  1=num spikes, 2=amplitude
    plotSettings.sortDir = 'descend';
    plotSettings.normalize = 0;

    plotSettings.graphPos = [253 581 600 156];
    plotSettings.tracePos = [147 281 1205 396];


    plotSettings.tolimit =  20;              % set to 0 if you want to plot all bulk & axial signal plots.
    %  set to -1 if you want equal # of control and
    %  mutant plots, the latter is better for
    %  comparison as bulk signals will haveS identical
    %  y-axis scaling. To plot a specific number of
    %  plots, set tolimit to that number

    plotSettings.overlayplots = 0;
    plotSettings.controlname ='wildtype-control';

    plotSettings.numColumns = 1;

    if plotsettings.OAS == 0
        if plotSettings.normalize == 1
            plotSettings.traceylimit = [-.2 1];
            plotSettings.peakthreshold = .1;
        else
            plotSettings.traceylimit = [0 8000];
            plotSettings.peakthreshold = 1000;
        end
        plotSettings.trimExperimentLength =0;
    else
        plotSettings.traceylimit = [0 20];
        plotSettings.peakthreshold = 3;
        plotSettings.axylimit = [0 45];
        plotSettings.trimExperimentLength =1;

    end

    %% plot single trace
    plotSettings.singleSpike = 1;
    plotSettings.spikeWindow = 10;

end
end