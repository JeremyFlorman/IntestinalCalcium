function [mtdata, wtdata, settings] = processWormdata(wormdata,settings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<2

    settings = returnPlotSettings();
end

wormdata = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\itr-1\itr-1_mergedData.mat";

normalize = settings.normalize;

if ischar(wormdata) || isstring(wormdata)
    mtdata = load(wormdata);
    mtdata = mtdata.wormdata;
elseif isstruct(wormdata)
    mtdata = wormdata;
end

if isfield(mtdata, 'include')
    mtdata = mtdata(logical([mtdata.include]));
end


if isfield(mtdata, 'controlData')
    wtdata = mtdata(1).controlData;
    if isfield(wtdata, 'include')
        wtdata = wtdata(logical([wtdata.include]));
    end
else
    wtdata = mtdata;
end


wtdata = subtractBackground(wtdata, settings);
mtdata = subtractBackground(mtdata, settings);



% normalize bulk signal by dividing by the mean wildtype bulk signal.
if normalize == 1
    wtbulksig = cell2mat({wtdata(:).bulkSignal});
    wtmedian = median(wtbulksig,'all','omitnan');
    for i = 1:length(mtdata)
        mtdata(i).bulkSignal = mtdata(i).bulkSignal./wtmedian-1;
    end

    for i = 1:length(wtdata)
        wtdata(i).bulkSignal = wtdata(i).bulkSignal./wtmedian-1;
    end
end




mtdata = processSpikes(mtdata,settings);
wtdata = processSpikes(wtdata,settings);


end
%% background subtract bulk and axial Signal
function [processedData] = subtractBackground(inputData, settings)
for i = 1:length(inputData)
    
    inputData(i).bulkSignal = fillmissing(inputData(i).bulkSignal-inputData(i).backgroundSignal, 'movmedian',100);

    backgroundMatrix = repmat(inputData(i).backgroundSignal,1,size(inputData(i).autoAxialSignal,2));
     axsig = inputData(i).autoAxialSignal-backgroundMatrix;

    if settings.autoFixAxialSignal
        toQuerry = settings.axSigToQuerry;
        inputData(i).autoAxialSignal = autoFixSignal(axsig,toQuerry);
    else 
        inputData(i).autoAxialSignal = axsig;
    end


end
processedData = inputData;
end




%% Spike processing

function [processedData] = processSpikes(inputData, settings)

peakdistance = settings.peakdistance;
peakthreshold = settings.peakthreshold;
peakwidth = settings.peakwidth;
sortType = settings.sortType;
sortDir = settings.sortDir;
secondsPrePost = settings.spikeProfileWindow;
framerate = settings.framerate;

timePreSpike = framerate*secondsPrePost;
timePostSpike = framerate*secondsPrePost;


num = nan(length(inputData),1);
AvAmp = nan(length(inputData),1);
sortOrder = nan(length(inputData),1);

if settings.showFitParams == 1
    figure()
    t = tiledlayout(ceil(length(inputData)/2),2,Padding='tight',TileSpacing='tight');
end

plotPropagationRate = 1;

if plotPropagationRate == 1
    figure();
    propAx = axes();
end


for i = 1:length(inputData)
    peakProfiles = [];
    tempint = [];
    bulkSignal = inputData(i).bulkSignal;
    axialSignal = inputData(i).autoAxialSignal;
    chunksize = floor(0.33*size(axialSignal,2));

    [tempamp, templocs] = findpeaks(bulkSignal, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance,'MinPeakWidth',peakwidth);

    if length(templocs)>1
        tempint = diff(templocs)/framerate;
    end

    if ~isempty(templocs)
        baseline = nan(length(templocs)+1,1);
        splitpoints = nan(length(templocs)+1,1);
        splitpoints(end) = length(bulkSignal);
        AUC = nan(length(templocs),1);
        propagationRate = nan(length(templocs),1);



        for q = 1:length(templocs)

            pre = templocs(q)-timePreSpike;
            post = templocs(q)+timePostSpike;

            if pre>0 && post<= length(bulkSignal) % spike traces
                ttrace = bulkSignal(pre:post);
                peakProfiles = horzcat(peakProfiles, ttrace);
%                 axialPeak = axialSignal(pre:post,:)';
            end


            %% Peak kinetics 

            if q==1 % find timepoints for valleys between peaks.
                segX = 1:templocs(q); % first valley
            else
                segX = templocs(q-1):templocs(q); % all other valleys
            end

            seg = bulkSignal(segX); % find values for valleys between peaks.
            segsort = sort(seg); % sort valley points from lowest to highest.
            baseline(q) = mean(segsort(1:floor(length(segsort)*.05))); % take the mean of the lowest 50% of values in the valley


            if q == 1
                splitpoints(q) = 1;
            else
                splitpoints(q) = segX(islocalmin(seg,'MaxNumExtrema',1));
            end

            % get the last baseline
            if q==length(templocs)
                lastSegX = templocs(q):length(bulkSignal);
                seg = bulkSignal(lastSegX); % find values for valleys between peaks.
                segsort = sort(seg); % sort valley points from lowest to highest.
                baseline(q+1) = mean(segsort(1:floor(length(segsort)*.05))); % take the mean of the lowest 50% of values in the valley
                if settings.showFitParams == 1
                line(lastSegX, repmat(baseline(q+1),length(lastSegX),1), 'Color', 'k', 'LineStyle', ':', 'Parent', fitAx)
                end
            end


            if settings.showFitParams == 1
                if q == 1
                    fitAx = nexttile(t);
                    plot(bulkSignal,'Parent', fitAx)
                    fitAx.XAxis.Visible = 0;
                    fitAx.YAxis.Visible = 0;
                    xlim([0 length(bulkSignal)])

                    if iscell(inputData(i).filename)
                       fn= inputData(i).filename{:};
                    elseif ischar(inputData(i).filename)
                        fn = inputData(i).filename;
                    end

                    expname = strsplit(fn,'\');
                    title(fitAx, expname{end},'FontSize',8,'Interpreter','none')
                end
                line(segX, repmat(baseline(q),length(segX),1), 'Color', 'k', 'LineStyle', '-','Parent', fitAx)
                line([splitpoints(q) splitpoints(q)], [0 1000],'Color', 'k', 'LineStyle', ':','Parent', fitAx )
            end
        end

        riseX = cell(length(templocs),1);
        riseY = cell(length(templocs),1);
        fallX = cell(length(templocs),1);
        fallY = cell(length(templocs),1);
        rTime = nan(length(templocs),1);
        fTime = nan(length(templocs),1);

        for j = 1:length(templocs) % split signal into rise and fall segments

            rX = splitpoints(j):templocs(j); % get X points for rise segment
            rY = bulkSignal(rX);  % get corresponding values for rise;

            fX = templocs(j):splitpoints(j+1); % get X points for fall segment
            fY = bulkSignal(fX);  % get corresponding values for fall;

            riseDiff = tempamp(j)- baseline(j); % find difference between rise baseline and peak.
            rise10 = baseline(j)+(riseDiff*0.1); % find 10% above baseline
            rise90 = tempamp(j)*0.9; % find 90% of max amplitude

            fallDiff = tempamp(j)-baseline(j+1); % find difference between peak and the fall baseline
            fall10 = baseline(j+1)+(fallDiff*0.1); % find 10% above baseline
            fall90 = tempamp(j)*0.9; % find 90% of max amplitude


            rstart = find(rY>rise10,1); % find first point
            rend = find(rY> rise90,1);

            if rstart>1
                riseX(j) = {rX(rstart:rend)};
                riseY(j) = {rY(rstart:rend)};
                rTime(j) = length(riseX{j})/settings.framerate;
            end



            fstart = find(fY<fall90,1);
            fend = find(fY<fall10, 1);
            fallX(j) = {fX(fstart:fend)};
            fallY(j) = {fY(fstart:fend)};
            fTime(j) = length(fallX{j})/settings.framerate;

            if settings.showFitParams == 1
                line(riseX{j},riseY{j}, 'Color', 'g', 'LineStyle', ':','Marker','^', 'MarkerSize', 2,'Parent', fitAx)
                line(fallX{j},fallY{j}, 'Color', 'r', 'LineStyle', ':','Marker', 'v','MarkerSize', 2,'Parent', fitAx)
            end

            %% wave propagation rate
            
            padValue = 45;          
            if rstart == 1 || riseX{j}(1)-padValue <1
                prePad = 1;
            else
                prePad = riseX{j}(1)-padValue;
            end

            postPad = riseX{j}(end)+padValue;
            if postPad > length(axialSignal)
                postPad = length(axialSignal);
            end
            
            axialPeak = axialSignal(prePad:postPad,:)'; % get axial signal during rising time

            head = mean(axialPeak(10:chunksize+10,:),1); % axial signal in head segment
            tail = mean(axialPeak(end-chunksize:end,:),1);  % axial signal in tail segment

            smWindow = 60;      % smoothing window/polynomial degree (sgolay method)
            smMethod = 'movmean';   % smoothing method
            smHead = smoothdata(head, smMethod,smWindow);
            smTail = smoothdata(tail, smMethod,smWindow);

            [~, headmax] = max(diff(smHead));
            [~, tailmax] = max(diff(smTail));

            propagationRate(q) = (headmax-tailmax)/settings.framerate;



            if plotPropagationRate == 1
                x = 1:size(head,2);
                plot(x,smHead,x,smTail,'Parent',propAx);
                title(['Propagation time: ' num2str(propagationRate(q)) ' seconds'])
                line([tailmax tailmax], propAx.YLim, 'linestyle', ':','Color', [.3 .3 .3], 'linewidth', 1)
                line([headmax headmax], propAx.YLim,'linestyle', ':', 'Color', [.3 .3 .3], 'linewidth' ,1)
                pause(1)
            end
%             imagesc(axialPeak, 'Parent', propAx)

        end


    



    elseif isempty(templocs)
        rTime = nan;
        fTime = nan;
        AUC = nan;
        propagationRate = nan;
    end



    riseNan = isnan(rTime);
    fallNan = isnan(fTime);
    aucNan = isnan(AUC); 

    inputData(i).peakTraces = peakProfiles;
    inputData(i).riseTime = rTime(~riseNan);
    inputData(i).fallTime = fTime(~fallNan);
    inputData(i).AUC = AUC(~aucNan);
    inputData(i).peakIntervals = tempint;
    inputData(i).peakAmplitude = tempamp;
    inputData(i).peakLoc = templocs;
%     inputData(i).propagationRate = propagationRate;

    if ~isempty(templocs)
        num(i) = length(templocs);
        AvAmp(i) = mean(tempamp);
    else
        num(i) = 0;
        AvAmp(i) = 0;
    end
end

if sortType == 0                % dont sort
    sortOrder = 1:length(inputData);
    if strcmpi(sortDir, 'ascend')
        flipud(sortOrder)
    end
elseif sortType == 1            % sort by number of spikes
    mtsorttype = num;
    [~, sortOrder] = sort(mtsorttype,sortDir);
elseif sortType == 2            % sort by spike amplitude
    mtsorttype = AvAmp;
    [~, sortOrder] = sort(mtsorttype,sortDir);
end

for i = 1:length(sortOrder)
    inputData(i).sortOrder = sortOrder(i);
end



processedData = inputData(sortOrder);

processedData(1).riseVector = vertcat(inputData(:).riseTime);
processedData(1).fallVector = vertcat(inputData(:).fallTime);
processedData(1).AUCVector = vertcat(inputData(:).AUC);
processedData(1).intervalVector = vertcat(inputData(:).peakIntervals);
% processedData(1).propagationVector = vertcat(inputData(:).propagationRate);
end





