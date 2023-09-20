function [] = plotAxialSignal(data,settings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

axylimit = settings.axylimit;
XTickInt = settings.axialXticint;
plotlimit = settings.tolimit;

if plotlimit == 0 || plotlimit>length(data) 
    num2plot = length(data);
else
    num2plot = plotlimit;
end


buffer = NaN(length(data(1).autoAxialSignal), 1);

for i = 1:num2plot
axsig = data(i).autoAxialSignal;
    if i == 1
        axialMatrix = axsig;
    elseif i>1
        axialMatrix = horzcat(axialMatrix,buffer, axsig);
    end
end


imagesc(smoothdata(axialMatrix,'movmedian',60)',axylimit)
colormap(viridis)

if isfield(data, 'genotype')
    title(['\it' data(1).genotype])
end

durationInMin = length(data(1).autoAxialSignal)/900;
framesInMin = length(data(1).autoAxialSignal)/durationInMin;

xt = [1 framesInMin*XTickInt:framesInMin*XTickInt:length(data(1).autoAxialSignal)];
xtl = 0:2:durationInMin;


ax = gca;
ax.XTick = xt;
ax.XTickLabel = xtl;
ax.YTickLabel = [];

xlabel('Time (min)');
end



