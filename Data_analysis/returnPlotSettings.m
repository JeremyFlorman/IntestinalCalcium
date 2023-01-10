function [plotSettings] = returnPlotSettings()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 %default 1000
plotSettings.peakdistance = 150;
plotSettings.peakwidth = 15;

plotSettings.axylimit = [0 25000]; % w/ background subtraction
% plotSettings.axylimit = [2000 35000]; % w/o background subtraction
plotSettings.wtcolor = [0.4 0.4 0.4];
plotSettings.mtcolor = [0.66 0.74 0.91];
plotSettings.mtedgecolor = [0.66 0.74 0.91] -.4;
plotSettings.loclinewidth = 1.5;
plotSettings.binedges = 0:2:150;     % for histograms/correlation
plotSettings.xlimits = [0 10];  % for traces 


%% sends input to 
% plotSettings.BulkSignalType =0; % 0 = raw bulkSignal 
%                                 % 1 = background subtracted signal (will
%                                 % throw error if "backgroundSignal field is
%                                 % absent.
%                                 % 3 = normalized to mean control signal

plotSettings.sortType = 1; % 1=num spikes, 2=amplitude
plotSettings.sortDir = 'descend';
plotSettings.normalize = 0;

switch plotSettings.normalize
    case 1
        plotSettings.traceylimit = [-.3 1.2];
        plotSettings.peakthreshold = .05;
    case 0 
        plotSettings.traceylimit = [-500 9000]; %w/ background
%         plotSettings.traceylimit = [4000 12000]; %no background 
        plotSettings.peakthreshold = 750;
end

end