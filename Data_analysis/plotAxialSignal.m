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


buffer = NaN(1,length(data(1).autoAxialSignal))';
stimY = NaN(num2plot,1);

for i = 1:num2plot
axsig = data(i).autoAxialSignal;
    if i == 1
        axialMatrix = axsig;
        stimY(i) = 5;
%         sizeIncrease = size(axialMatrix,2);
    elseif i>1
        axialMatrix = horzcat(axialMatrix,buffer, axsig);
        stimY(i) = stimY(i-1)+size(axsig,2)+1;
    end 
end


imagesc(smoothdata(axialMatrix,'movmedian',60)',axylimit)
colormap("turbo")

if isfield(data, 'stimTimes')
    hold on
    for i =1:num2plot
        stimtimes = data(i).stimTimes; 
        yvals = repmat(stimY(i), length(stimtimes),1);
        plot(stimtimes,yvals,'v','color' ,[0 0 0],'MarkerFaceColor',[.8 .8 .8], 'MarkerSize',8)
    end
end
hold off

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
box off
xlabel('Time (min)');
end



