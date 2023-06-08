function [] = plotOverlay(data, settings)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

yyaxis left
plotAxialSignal(data, settings)

yyaxis right


overlayBulkSignal(data, settings)
ax = gca;

ax.YColor = [1 1 1];
box off

end