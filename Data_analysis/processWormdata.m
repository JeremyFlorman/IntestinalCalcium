function [mtdata, wtdata, settings] = processWormdata(wormdata,settings)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<2

    settings = returnPlotSettings();
end

% wormdata = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\dec-9\dec-9_mergedData.mat";

normalize = settings.normalize;
trimExperimentLength = settings.trimExperimentLength;



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

if trimExperimentLength == 1
    mtdata = trimExp(mtdata);
    wtdata = trimExp(wtdata);
end


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
validatePropagationRate = settings.validatePropagationRate;
validateRiseFall = settings.validateRiseFall;
propMethod =1;


timePreSpike = framerate*secondsPrePost;
timePostSpike = framerate*secondsPrePost;


num = nan(length(inputData),1);
AvAmp = nan(length(inputData),1);
sortOrder = nan(length(inputData),1);

if settings.showFitParams == 1
    figure()
    t = tiledlayout(ceil(length(inputData)/2),2,Padding='tight',TileSpacing='tight');
end


if validatePropagationRate == 1
    figure();
    axFig = tiledlayout(2,1);
    propAx = nexttile();
    axAx = nexttile();
end

if validateRiseFall == 1
    figure();
    rAx = axes;
    keepValidating = 1;
end


for i = 1:length(inputData)
    peakProfiles = [];
    tempint = [];
    bulkSignal = inputData(i).bulkSignal;
    axialSignal = inputData(i).autoAxialSignal;
    chunksize = floor(0.2*size(axialSignal,2));

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


            %% wave propagation rate
            axPre = templocs(q)-floor(timePreSpike/3);
            axPost = templocs(q)+floor(timePreSpike/3);
            if axPre>0 && axPost<= length(axialSignal)
                try
                    axialPeak = smoothdata(axialSignal(axPre:axPost,:)',2, 'movmean', 60,'omitnan');

                    headstart = 25;
                    head = mean(axialPeak(headstart:chunksize+headstart,:),1); % axial signal in head segment
                    tail = mean(axialPeak(end-chunksize:end,:),1);  % axial signal in tail segment
                    [hpk, hloc] = findpeaks(head, 'SortStr','descend');
                    [tpk, tloc] = findpeaks(tail, 'SortStr', 'descend');

                    headRise = head(1:hloc(1))';
                    tailRise = tail(1:tloc(1))';

                    % find inflection point using derivative
                    [~, headmax] = max(diff(headRise));
                    [~, tailmax] = max(diff(tailRise));

                    % find inflection point using full-width@half-maximum
                    hFWHM = find(headRise<= hpk(1)/2,1,'last');
                    tFWHM = find(tailRise<= tpk(1)/2,1,'last');


                    if propMethod == 1 % derivative
                        propagationRate(q) = (headmax-tailmax)/settings.framerate;
                        hInflect = headmax;
                        tInflect = tailmax;

                    elseif propMethod == 2 % full width @ half maximum
                        hInflect = hFWHM;
                        tInflect = tFWHM;
                        if ~isempty(hFWHM) && ~isempty(tFWHM)
                            propagationRate(q) = (hFWHM-tFWHM)/settings.framerate;
                        end
                    elseif propMethod == 3 % peak location
                        hInflect = hloc(1);
                        tInflect = tloc(1);
                        if ~isempty(hloc) && ~isempty(tloc)
                            propagationRate(q) = (hloc(1)-tloc(1))/settings.framerate;
                        end
                    end

                catch
                end


                    if validatePropagationRate == 1
                        if propagationRate(q)>5 || propagationRate(q) <0
                            x = 1:size(head,2);
                            plot(x,head,x,tail,'Parent',propAx);

                            hold(propAx, "on")

                            plot(hloc(1),hpk(1)*1.05, 'v', 'MarkerFaceColor', [.07 .62 1],'MarkerEdgeColor', [.07 .62 1],'Parent',propAx)
                            if ~isempty(hInflect)
                                plot(hInflect,head(hInflect)*1.05, 'v', 'MarkerFaceColor', [.07 .62 1],'MarkerEdgeColor', [.07 .62 1],'Parent',propAx)
                            end

                            plot(tloc(1),tpk(1)*1.05, 'v', 'MarkerFaceColor', [.93 .69 .13],'MarkerEdgeColor', [.93 .69 .13],'Parent',propAx)
                            if ~isempty(tInflect)
                                plot(tInflect,tail(tInflect)*1.05, 'v', 'MarkerFaceColor', [.93 .69 .13],'MarkerEdgeColor', [.93 .69 .13],'Parent',propAx)
                            end

                            hold(propAx, "off")
                            xlim(propAx, [0 length(head)])
                            %
                            %                     line([headmax headmax], propAx.YLim,'linestyle', ':', 'Color', [.07 .62 1], 'linewidth' ,1.5, 'Parent', propAx)
                            %                     line([tailmax tailmax], propAx.YLim, 'linestyle', ':','Color', [.93 .69 .13], 'linewidth', 1.5, 'Parent', propAx)
                            %
                            imagesc(axialPeak,'Parent', axAx)
                            pkSz = size(axialPeak);
                            rectangle('Position',[0 pkSz(1)-chunksize pkSz(2) chunksize],'linestyle', ':','EdgeColor', [.93 .69 .13], 'linewidth', 1.5, 'Parent', axAx)
                            rectangle('Position',[0 headstart pkSz(2), chunksize],'linestyle', ':','EdgeColor', [.07 .62 1], 'linewidth', 1.5, 'Parent', axAx)

                            line([hInflect hInflect], [headstart headstart+chunksize], 'Color', [.07 .62 1], 'linewidth' ,1.5, 'Parent', axAx)
                            line([tInflect tInflect], [pkSz(1) pkSz(1)-chunksize],'Color', [.93 .69 .13], 'linewidth', 1.5, 'Parent', axAx)

                            title(['Propagation time: ' num2str(propagationRate(q)) ' seconds'], 'Parent',axFig)



                            txt = input("Look ok? if not press letter key before hitting enter. Type 'exit' to quit","s");

                            if ~isempty(txt)
                                propagationRate(q) = NaN;

                            end

                            if ~isempty(txt) && strcmp(txt,'exit')
                                validatePropagationRate = 0;
                            end

                        end

                    end



            end

            %% Peak kinetics

            pre = templocs(q)-timePreSpike;
            post = templocs(q)+timePostSpike;

            if pre>0 && post<= length(bulkSignal) % spike traces
                ttrace = bulkSignal(pre:post);
                peakProfiles = horzcat(peakProfiles, ttrace);
            end

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
                    line(lastSegX, repmat(baseline(q+1),length(lastSegX),1), 'Color', 'k', 'LineStyle', ':')
                end
            end


            if settings.showFitParams == 1
                if q == 1
                    nexttile(t)
                    plot(bulkSignal)
                    ax = gca;
                    ax.XAxis.Visible = 0;
                    ax.YAxis.Visible = 0;
                    xlim([0 length(bulkSignal)])

                    if iscell(inputData(i).filename)
                        fn= inputData(i).filename{:};
                    elseif ischar(inputData(i).filename)
                        fn = inputData(i).filename;
                    end

                    expname = strsplit(fn,'\');
                    title(ax, expname{end},'FontSize',8,'Interpreter','none')
                end
                line(segX, repmat(baseline(q),length(segX),1), 'Color', 'k', 'LineStyle', '-')
                line([splitpoints(q) splitpoints(q)], [0 1000],'Color', 'k', 'LineStyle', ':' )
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
                line(riseX{j},riseY{j}, 'Color', 'g', 'LineStyle', ':','Marker','^', 'MarkerSize', 2)
                line(fallX{j},fallY{j}, 'Color', 'r', 'LineStyle', ':','Marker', 'v','MarkerSize', 2)
            end

            if settings.validateRiseFall == 1
                if keepValidating == 1
                    if j == 1
                        plot(bulkSignal, Parent=rAx)
                    end

                    line(riseX{j},riseY{j}, 'Color', 'g', 'LineStyle', ':','Marker','^', 'MarkerSize', 2)
                    if ~isempty(riseX{j})
                        rtxt = input('Rise OK? (y/n)... type exit to quit','s');
                        if strcmpi(rtxt, 'n')
                            rTime(j) = nan;
                        elseif strcmpi(rtxt, 'exit')
                            keepValidating = 0;
                            break
                        end
                    end

                    line(fallX{j},fallY{j}, 'Color', 'r', 'LineStyle', ':','Marker', 'v','MarkerSize', 2)
                    if ~isempty(fallX{j})
                        ftxt = input('fall OK? (y/n)... type exit to quit','s');
                        if strcmpi(ftxt, 'n')
                            fTime(j) = nan;
                        elseif strcmpi(ftxt, 'exit')
                            keepValidating = 0;
                            break
                        end

                    end
                end
            end



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
    inputData(i).propagationRate = propagationRate;

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
        %         sortOrder = fliplr(sortOrder);
        sortOrder = sortOrder(randperm(length(sortOrder))); %shuffle 
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
processedData(1).propagationVector = vertcat(inputData(:).propagationRate);

if validateRiseFall == 1 && keepValidating == 1
    [file,path] = uiputfile('*.xlsx');
    riseFallSaveName = fullfile(path,file);
    writematrix(processedData(1).riseVector, riseFallSaveName, 'Sheet', 'Rise Time');
    writematrix(processedData(1).fallVector, riseFallSaveName, 'Sheet', 'Fall Time');
end

end

%% trim experiments
function [processedData] = trimExp(inputData)
lens = nan(length(inputData),1);

for i = 1:length(inputData)
    lens(i)=length(inputData(i).bulkSignal);
end

minlen = min(lens);

for i = 1:length(inputData)
    inputData(i).autoAxialSignal = inputData(i).autoAxialSignal(1:minlen,:);
    inputData(i).sumSignal = inputData(i).sumSignal(1:minlen);
    inputData(i).bulkSignal = inputData(i).bulkSignal(1:minlen);
    inputData(i).bulkAboveBkg = inputData(i).bulkAboveBkg(1:minlen);
    inputData(i).backgroundSignal = inputData(i).backgroundSignal(1:minlen);
    inputData(i).orientation = inputData(i).orientation(1:minlen);
    inputData(i).area = inputData(i).area(1:minlen);
    if isfield(inputData,'velocity')
        inputData(i).velocity = inputData(i).velocity(1:minlen);
    end

    if isfield(inputData,'stimTimes')
        validtimes = inputData(i).stimTimes<minlen;
        inputData(i).stimTimes = inputData(i).stimTimes(validtimes);
    end
    processedData = inputData;

end


end



