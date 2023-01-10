function [] = plotOverlay(datapath, plotcontrol,plotlimit)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

yyaxis left
plotAxialSignal(datapath,plotcontrol,plotlimit)

yyaxis right


overlayBulkSignal(datapath, plotcontrol,plotlimit)
ax = gca;

ax.YColor = [1 1 1];
box off

end