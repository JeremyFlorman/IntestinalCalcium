function [] = plotAxialSignal(data,settings,labelXAxis)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

axylimit = settings.axylimit;
XTickInt = settings.axialXticint;
plotlimit = settings.tolimit;
axSigCMap = settings.axSigCMap;
framerate = settings.framerate;

if plotlimit == 0 || plotlimit>length(data) 
    num2plot = length(data);
else
    num2plot = plotlimit;
end


% buffer = NaN(1,length(data(1).autoAxialSignal))';
buffer = repmat(axylimit(1)-1, length(data(1).autoAxialSignal),1);
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
% imagesc(imgaussfilt(axialMatrix', 2), axylimit)
colormap(axSigCMap);

if isfield(data, 'stimTimes')
    hold on
    for i =1:num2plot
        stimtimes = data(i).stimTimes; 
        yvals = repmat(stimY(i), length(stimtimes),1);
        plot(stimtimes,yvals,'v','color' ,[0 0 0],'MarkerFaceColor',[.8 .8 .8], 'MarkerSize',4)
    end
end
hold off

%% Annotate pBocs
if isfield(data, 'pBoc')
    for i = 1:num2plot
        bocTimes = data(i).pBoc;
        bocY = stimY(i);
        for k = 1:length(bocTimes)
            text(bocTimes(k), bocY, 'pBoc', 'Color',[1 0.9 0.5], 'FontSize', 5)
        end
    end
end


if settings.annotateFood == 1
    for i =1:num2plot
        if isfield(data(i), 'onFood')
            for k = 1:length(data(i).onFood)
                foodStart = data(i).onFood(k);

                if k<=length(data(i).offFood)
                    foodEnd = data(i).offFood(k);
                else
                    foodEnd = find(~isnan(data(i).bulkSignal),1,'last');
                end
                foodY = stimY(i)+5;
                line([foodStart foodEnd], [foodY foodY], 'Color', [0.97 0.93 0.62], 'LineWidth', 0.5, 'LineStyle', '-', 'Marker', 'none')
            end
        end
    end
end


if isfield(data, 'genotype')
    title(['\it' data(1).genotype])
end

durationInMin = length(data(1).autoAxialSignal)/(framerate*60);
framesInMin = length(data(1).autoAxialSignal)/durationInMin;

xt = [1 framesInMin*XTickInt:framesInMin*XTickInt:length(data(1).autoAxialSignal)];
xtl = 0:XTickInt:durationInMin;


ax = gca;
ax.XTick = xt;
ax.XTickLabel = xtl;
ax.YTickLabel = [];
% ax.YAxis.Visible = 'off'
box off

if labelXAxis ==1
xlabel('Time (min)')
else
    ax.XTickLabel = [];
end

end