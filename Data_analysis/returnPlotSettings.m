function [plotSettings] = returnPlotSettings()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 %default 1000
plotSettings.peakdistance = 150;
plotSettings.peakwidth = 15;
plotSettings.autoFixAxialSignal = 1;
plotSettings.axSigToQuerry = 0.15;
plotSettings.framerate = 15;

plotSettings.axylimit = [-500 12000];
plotSettings.wtcolor = [0.4 0.4 0.4];
plotSettings.mtcolor = [0.66 0.74 0.91];
plotSettings.mtedgecolor = [0.66 0.74 0.91] -.4;
plotSettings.loclinewidth = 1.5;
plotSettings.binedges = 0:2:90;     % for histograms/correlation
plotSettings.xlimits = [0 10];  % for traces 
plotSettings.spikeProfileWindow = 30;

plotSettings.axialXticint = 2; % # x tick interval (in minutes) for axial signal plots. 
%% sends input to 
% plotSettings.BulkSignalType =0; % 0 = raw bulkSignal 
%                                 % 1 = background subtracted signal (will
%                                 % throw error if "backgroundSignal field is
%                                 % absent.
%                                 % 3 = normalized to mean control signal

plotSettings.sortType = 0; % 0=dont sort,  1=num spikes, 2=amplitude
plotSettings.sortDir = 'descend';
plotSettings.normalize = 0;

switch plotSettings.normalize
    case 1
        plotSettings.traceylimit = [-.2 1];
        plotSettings.peakthreshold = .1;
    case 0 
        plotSettings.traceylimit = [-500 5000]; % Rebekka 220112
%         plotSettings.traceylimit = [-500 9000]; %w/ background
%         plotSettings.traceylimit = [4000 12000]; %no background 
        plotSettings.peakthreshold = 500;
end

end