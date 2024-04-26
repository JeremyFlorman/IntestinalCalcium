function [] = plotOverlay(data, settings,labelXAxis)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

yyaxis left
plotAxialSignal(data, settings,labelXAxis)
ax = gca;
ax.YAxis(1).Visible = "off";

yyaxis right


overlayBulkSignal(data, settings,labelXAxis)
ax = gca;

ax.YAxis(2).Visible = "off";
box off

end