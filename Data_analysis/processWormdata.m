function [mtdata, wtdata, settings] = processWormdata(wormdata,settings)
%PROCESSWORMDATA Process worm data for analysis
%   [MTDATA, WTDATA, SETTINGS] = PROCESSWORMDATA(WORMDATA, SETTINGS)
%   processes the worm data provided in WORMDATA using the settings
%   specified in SETTINGS. If SETTINGS is not provided, default settings
%   are used. The function returns the processed mutant data (MTDATA),
%   wild-type data (WTDATA), and the settings used (SETTINGS).
%
%   Inputs:
%       wormdata - A file path (string or character array) to the data or a
%                  structure containing the data.
%       settings - A structure containing various settings for data
%                  processing. If not provided, default settings are used.
%
%   Outputs:
%       mtdata   - Processed mutant data.
%       wtdata   - Processed wild-type data.
%       settings - The settings used for processing.
%
%   Description:
%       1. The function starts by checking if SETTINGS are provided. If not,
%          it uses default settings from RETURNPLOTSETTINGS.
%
%       2. If WORMDATA is a file path (string or character array), it loads
%          the data from the specified file. If WORMDATA is already a
%          structure, it uses it directly.
%
%       3. It extracts control data (WTDATA) from the mutant data (MTDATA)
%          if the CONTROLDATA field is present.
%
%       4. The function filters the data to include only the relevant entries
%          based on the INCLUDE field.
%
%       5. The SUBTRACTBACKGROUND function is called to remove background
%          signals from the bulk and axial signals in the data. This function
%          handles normalization of axial signals if specified.
%
%       6. Depending on the NORMALIZE setting in SETTINGS, the function
%          normalizes the data:
%           - If NORMALIZE is 2, it normalizes the data using z-scores.
%           - If NORMALIZE is 3, it normalizes the bulk signal by dividing by
%             the mean wild-type bulk signal.
%
%       7. If specified in the settings (TRIMEXPERIMENTLENGTH or
%          ANALYZEPARTIAL), the function trims the experiment length to focus
%          on specific parts of the data.
%
%       8. The PROCESSSPIKES function is called to detect and process spikes
%          in the bulk and axial signals. It calculates various parameters
%          such as peak distance, peak threshold, and propagation rates. It
%          also validates spike characteristics if specified in the settings.
%
%   Example:
%       settings = returnPlotSettings();
%       [mtdata, wtdata, settings] = processWormdata('path/to/data.mat', settings);
%
%   See also RETURNPLOTSETTINGS
%   Copyright 2024 Jeremy Florman

if nargin<1
    settings = returnPlotSettings();
    wormdata = evalin("caller",'wormdata');
end

% wormdata = "C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\dec-9\dec-9_mergedData.mat";

normalize = settings.normalize;
trimExperimentLength = settings.trimExperimentLength;
analyzePartial = settings.analyzePartial;
saveWormdata2workspace = settings.saveWormdata2workspace;
trimstim = settings.trim2stim;
trimSpike = settings.trim2spike;
trimOnFood = settings.trim2onFood;
trimOffFood = settings.trim2offFood;


if ischar(wormdata) || isstring(wormdata)
    mtdata = load(wormdata);
    mtdata = mtdata.wormdata;
elseif isstruct(wormdata)
    mtdata = wormdata;
end


if isfield(mtdata, 'controlData')
    for i = 1:length(mtdata) % find a non-empty cell
        if~isempty(mtdata(i).controlData)
            wtdata = mtdata(i).controlData;
        end
    end
    if isfield(wtdata, 'include')
        wtdata = wtdata(logical([wtdata.include]));
    end
else
    if isfield(mtdata, 'include')
        mtdata = mtdata(logical([mtdata.include]));
    end
    wtdata = mtdata;
end

if isfield(mtdata, 'include')
    mtdata = mtdata(logical([mtdata.include]));
end

wtdata = subtractBackground(wtdata, settings);
mtdata = subtractBackground(mtdata, settings);


if strcmp(normalize, 'Delta F/F0') == 1 % normalize by deltaF/F0
    mtdata = deltaF(mtdata,settings);
    wtdata = deltaF(wtdata,settings);
elseif strcmp(normalize, 'Z-Score') == 1 % normalize by calculating z-score
    for i = 1:length(mtdata)
        mtdata(i).bulkSignal = zscore(mtdata(i).bulkSignal);
    end

    for i = 1:length(wtdata)
        wtdata(i).bulkSignal = zscore(wtdata(i).bulkSignal);
    end

elseif strcmp(normalize, 'Control') == 1 % normalize bulk signal by dividing by the mean wildtype bulk signal.
    wtbulksig = cell2mat({wtdata(:).bulkSignal});
    wtmedian = median(wtbulksig,'all','omitnan');
    for i = 1:length(mtdata)
        mtdata(i).bulkSignal = mtdata(i).bulkSignal./wtmedian-1;
    end

    for i = 1:length(wtdata)
        wtdata(i).bulkSignal = wtdata(i).bulkSignal./wtmedian-1;
    end
end




%% Align experiments to event
if trimstim || trimSpike || trimOnFood || trimOffFood 
    wtdata = alignTraces(wtdata,settings);
    mtdata = alignTraces(mtdata,settings);
end


if trimExperimentLength == 1 || analyzePartial == 1
    mtdata = trimExp(mtdata,settings);
    wtdata = trimExp(wtdata, settings);
else
    mtdata = padTraces(mtdata, settings);
    wtdata = padTraces(wtdata, settings);
end


mtdata = processFood(mtdata, settings);
wtdata = processFood(wtdata, settings);

mtdata = processSpikes(mtdata,settings, 0);
wtdata = processSpikes(wtdata,settings, 1);


%%  Concatenating fields
mtdata(1).int1Matrix = horzcat(mtdata(:).int1Signal);
mtdata(1).int9Matrix = horzcat(mtdata(:).int9Signal);

wtdata(1).int9Matrix = horzcat(wtdata(:).int9Signal);
wtdata(1).int9Matrix = horzcat(wtdata(:).int9Signal);

if isfield(mtdata, 'pumpingRate')
    mtdata(1).pumpingMatrix = horzcat(mtdata(:).pumpingRate);
    wtdata(1).pumpingMatrix = horzcat(wtdata(:).pumpingRate);
end

mtdata(1).timeInSec = linspace(0, length(mtdata(1).bulkSignal)/settings.framerate,length(mtdata(1).bulkSignal))';
wtdata(1).timeInSec = linspace(0, length(wtdata(1).bulkSignal)/settings.framerate,length(wtdata(1).bulkSignal))';

%% loading processed data into base workspace
if saveWormdata2workspace == 1
    if isfield(mtdata, 'genotype')
        dataName = strrep(mtdata(1).genotype,'-','');
        dataName = strrep(dataName, '+', 'Plus');
        dataName = strrep(dataName, ' ', '');
        dataName = strrep(dataName, '(', '');
        dataName = strrep(dataName, ')', '');
        dataName = strrep(dataName, ';', '');

        assignin("base", [dataName 'Data'], mtdata)
        assignin("base", [dataName '_ControlData'], wtdata)
    else
        assignin("base", 'SingleSpikeData', mtdata);
    end
end



% assignin('base','normWT', wtdata(1).amplitudeVector/mean(wtdata(1).amplitudeVector,"omitmissing"))
% assignin('base', 'normMT', mtdata(1).amplitudeVector/mean(wtdata(1).amplitudeVector,"omitmissing"))

% for i = 1:length(mtdata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%% background subtract bulk and axial Signal
function [processedData] = subtractBackground(inputData, settings)
normAx = 0; % nomralize axial signal?

for i = 1:length(inputData)
    % Subtract Background and fill missing datapoints in bulk signal
    bkgSignal = smoothdata(inputData(i).backgroundSignal, 'movmedian',10*settings.framerate);

    if ~strcmp(settings.normalize, 'Delta F/F0')
        inputData(i).bulkSignal = smoothdata(fillmissing(inputData(i).bulkSignal-bkgSignal, 'movmedian',100), 'movmedian', 15);
    end

    % fill outliers more than 3 standard deviations outside moving mean
    % with previous non-outlier value
    % inputData(i).bulkSignal = filloutliers(inputData(i).bulkSignal,"spline","movmean",25,ThresholdFactor=3);

    % Subtract Background  in axial signal
    backgroundMatrix = repmat(bkgSignal,1,size(inputData(i).autoAxialSignal,2));
    axsig = inputData(i).autoAxialSignal-backgroundMatrix;

    if settings.autoFixAxialSignal && ~isfield(inputData(i), 'noAutoFix')
        toQuerry = settings.axSigToQuerry;
        inputData(i).autoAxialSignal = autoFixSignal(axsig,toQuerry);
    else
        inputData(i).autoAxialSignal = axsig;
    end

    if normAx == 1
        as = inputData(i).autoAxialSignal;
        axMean = smoothdata(median(as,1,'omitnan'),"gaussian", 20);
        inputData(i).autoAxialSignal = inputData(i).autoAxialSignal./axMean;
    end

end
processedData = inputData;
end

%% Spike processing
function [processedData] = processSpikes(inputData, settings, isControlData)


peakdistance = settings.peakdistance;
peakthreshold = settings.peakthreshold;
peakwidth = settings.peakwidth;
sortType = settings.sortType;
sortDir = settings.sortDir;
secondsPrePost = settings.spikeProfileWindow;
framerate = settings.framerate;
analyzePropagationRate = settings.analyzePropagationRate;
validatePropagationRate = settings.validatePropagationRate;
validateRiseFall = settings.validateRiseFall;
waveAverage = settings.waveAverage;

if isControlData == 1 
    analyzePropagationRate = 0; % don't run wave propagation analysis on control datasets
end


% pixel Scaling
if settings.isOAS == 1
    pxPerMM = 381; % pixel scaling for 4x objective in 2x2 binning on OAS
elseif settings.isOAS == 0
    pxPerMM = 200.071; % pixel scaling for 5x objective in 4x4 binning on Zeiss inverted scope
end

umPerPixel = 1000/pxPerMM;

timePreSpike = framerate*secondsPrePost;
timePostSpike = framerate*secondsPrePost;

num = nan(length(inputData),1);
AvAmp = nan(length(inputData),1);
sortOrder = nan(length(inputData),1);

if settings.showFitParams == 1
    figure()
    t = tiledlayout(ceil(length(inputData)/2),2,Padding='tight',TileSpacing='tight');
end


if analyzePropagationRate == 1 && validatePropagationRate == 1
    figure(Position=[888.2000 200.2000 435.2000 504.8000], Color=[1 1 1]);
end


% if validatePropagationRate == 1
%     figure('WindowStyle','docked');
%     axFig = tiledlayout(2,1);
%     propAx = nexttile();
%     axAx = nexttile();
% end

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

    if isfield(inputData, 'wormLength')
        wormLength = inputData(i).wormLength;
    end

    % chunksize = floor(0.2*size(axialSignal,2));

    [tempamp, templocs] = findpeaks(bulkSignal, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance,'MinPeakWidth',peakwidth);


    if length(templocs)>1
        tempint = diff(templocs)/framerate;
    end

    if ~isempty(templocs)
        baseline = nan(length(templocs)+1,1);
        splitpoints = nan(length(templocs)+1,1);
        splitpoints(end) = length(bulkSignal);
        axialPeak = nan(size(axialSignal,2), floor(timePreSpike/3)*2+1, length(templocs));

        if waveAverage == 0
        propagationRate = nan(length(templocs),settings.numSegments);
        rSquared = nan(length(templocs),settings.numSegments);
        validFlags = nan(length(templocs),settings.numSegments);
        end







        for q = 1:length(templocs)
            axPre = templocs(q)-floor(timePreSpike/3);
            axPost = templocs(q)+floor(timePreSpike/3);

            %% Wave Propagation Rate

            if axPre>0 && axPost<= length(axialSignal)
                axialPeak(:,:,q) = axialSignal(axPre:axPost,:)';

                if analyzePropagationRate == 1 &&  waveAverage == 0
                    filename = strrep(stripPath(inputData(i)), '_', ' ');
                    figTitle = [filename{:} ' spike # ' num2str(q)];
                    [propRate, r2, flags, fitImage] = getWavePropagationRate(axialPeak(:,:,q), mean(wormLength(axPre:axPost),'omitnan'), settings, figTitle);
                    propagationRate(q, 1:length(propRate)) = propRate;
                    rSquared(q, 1:length(r2)) = r2;
                    validFlags(q, 1:length(flags)) = flags;

                    if i == 1
                        fitMontage = fitImage;
                    else
                        fitMontage = horzcat(fitMontage, fitImage);
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
        tau = nan(length(templocs),1);

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

            %% Decay Constant
            isTau = 0;
            if ~isempty(fstart) && ~isempty(fend) && fend > fstart + 5
                isTau = 1;
                fallStart = fX(fstart);
                fallEnd   = fX(fend);

                t_fit = (0:(fallEnd - fallStart))' / settings.framerate; % time vector for decay curve
                y_fit = bulkSignal(fallStart:fallEnd); % signal values for decay curve

                % Fit exponential: A * exp(-t/tau) + C
                A0 = y_fit(1) - y_fit(end);
                tau0 = 2;  % seconds
                C0 = y_fit(end);
                b0 = [A0, tau0, C0];

                expModel = @(b, t) b(1) * exp(-t / b(2)) + b(3);


                opts = optimset('Display', 'off');
                bFit = fminsearch(@(b) sum((expModel(b, t_fit) - y_fit).^2), b0, opts);
                tau(j) = bFit(2);

                if tau(j)>1000
                    % figure();
                    % tau(j)
                    % plot(bulkSignal)
                    % hold on
                    % plot(fallStart:fallEnd, expModel(bFit, t_fit), 'm-', 'LineWidth', 1.2);
                    % hold off
                    % drawnow();
                    tau(j) = [];
                end


                if settings.showFitParams == 1
                    plot(fallStart:fallEnd, expModel(bFit, t_fit), 'm-', 'LineWidth', 1.2);
                end
            end


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

                    if isTau == 1
                        hold on
                        plot(fallStart:fallEnd, expModel(bFit, t_fit), 'm-', 'LineWidth', 1.2);
                        hold off
                    end

                    if ~isempty(fallX{j})
                        ftxt = input('fall OK? (y/n)... type exit to quit','s');
                        if strcmpi(ftxt, 'n')
                            fTime(j) = nan;
                            tau(j) = nan;
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
        % propagationRate = nan;
        axialPeak = nan;
        tau = nan;
    end
    
    %% Int1 and Int9 signal
    dataSize = size(inputData(1).autoAxialSignal);
    segmentSize = round(dataSize(2)*0.025);
    lengthChunk = floor(dataSize(2)/3);
    
    
    % Int1 start and end points
    anteriorMean = mean(inputData(i).autoAxialSignal(:,1:lengthChunk),1,'omitmissing');
    [~,antLoc] = findpeaks(rescale(anteriorMean), NPeaks=1, SortStr="descend");
    
    antStart = antLoc-segmentSize;
    antEnd = antLoc+segmentSize;
    if antStart<1
        antStart = 1;
    end

    % Int9 start and end points
    posteriorMean = mean(inputData(i).autoAxialSignal(:,end-lengthChunk:end),1,'omitmissing');
    [~,postLoc] = findpeaks(rescale(posteriorMean), NPeaks=1, SortStr="descend");

    postLoc = postLoc+(lengthChunk*2); % account for the anterior two thirds
    postStart = postLoc-segmentSize;
    postEnd = postLoc+segmentSize;
    if postEnd > dataSize(2)
        postEnd = dataSize(2);
    end

    %% plot int1/int9 identification
    % plot(mean(inputData(i).autoAxialSignal,1,'omitmissing'))
    % ax = gca;
    % lims = ax.YLim;
    % line([antLoc-segmentSize antLoc-segmentSize], [lims(1) lims(2)])
    % line([antLoc+segmentSize antLoc+segmentSize], [lims(1) lims(2)])
    % line([postLoc-segmentSize postLoc-segmentSize],  [lims(1) lims(2)])
    % line([postLoc+segmentSize postLoc+segmentSize],  [lims(1) lims(2)])
    %

    inputData(i).int1Signal = mean(inputData(i).autoAxialSignal(:,antStart:antEnd),2, 'omitmissing');
    inputData(i).int9Signal = mean(inputData(i).autoAxialSignal(:,postStart:postEnd),2, 'omitmissing');
    
    


    tau(tau==0) = nan;
    riseNan = isnan(rTime);
    fallNan = isnan(fTime);
    tauNan = isnan(tau);

    inputData(i).peakTraces = peakProfiles;
    inputData(i).riseTime = rTime(~riseNan);
    inputData(i).fallTime = fTime(~fallNan);
    inputData(i).tau = tau(~tauNan);
    inputData(i).peakIntervals = tempint;
    inputData(i).meanInterval = mean(tempint, 'omitmissing');
    inputData(i).peakAmplitude = tempamp;
    inputData(i).peakLoc = templocs;


    inputData(i).axialPeak = axialPeak;
    inputData(i).avgKymograph = mean(axialPeak,3,'omitmissing');

    %% calculate average propagation rate of spikes in a recording instead of individual spikes
    if analyzePropagationRate == 1 &&  waveAverage == 1
        if isfield(inputData, 'wormLength') && ~isempty(inputData(i).wormLength)
            filename = strrep(stripPath(inputData(i)), '_', ' ');
            figTitle = [filename{:} ' spike average'];
            [propagationRate, rSquared, validFlags, fitImage] = getWavePropagationRate(inputData(i).avgKymograph, mean(wormLength,'omitnan'), settings,figTitle);

            if i == 1
                fitMontage = fitImage;
            else
                fitMontage = horzcat(fitMontage, fitImage);
            end
       
        end
    end


    if analyzePropagationRate == 1
        inputData(i).propagationRate = propagationRate;
        inputData(i).rSquared = rSquared;
        inputData(i).validFlags = validFlags;
      
        if settings.isOAS == 1
            umPerPixel = 2.6247; % pixel scaling for 4x objective in 2x2 binning on OAS
        elseif settings.isOAS == 0
            umPerPixel = 1/199.641*1000; % pixel scaling for 5x objective in 4x4 binning on Zeiss inverted scope
        end

        inputData(i).medianWormLength = median(inputData(i).wormLength, 'omitmissing')*umPerPixel;
        inputData(i).medianAbsPropRate = median(abs(propagationRate), 'omitmissing');
    end

    if ~isempty(templocs)
        num(i) = length(templocs);
        AvAmp(i) = mean(tempamp);
        firstSpike(i) = templocs(1);
    else
        firstSpike(i) = 0;
        num(i) = 0;
        AvAmp(i) = 0;
    end

    meanSignal(i) = mean(inputData(i).bulkSignal,"all", "omitmissing");


end


%% sorting
if sortType == 0                % dont sort
    sortOrder = 1:length(inputData);
    if strcmpi(sortDir, 'ascend')
        sortOrder = fliplr(sortOrder);
    elseif strcmpi(sortDir, 'shuffle')
        sortOrder = sortOrder(randperm(length(sortOrder))); %shuffle
    end

elseif sortType == 1            % sort by number of spikes
    if strcmpi(sortDir, 'shuffle')
        disp('Shuffle only works with "dont sort", using descending order')
        sortDir = 'descend';
    end
    mtsorttype = num;
    [~, sortOrder] = sort(mtsorttype,sortDir);

elseif sortType == 2            % sort by spike amplitude
    if strcmpi(sortDir, 'shuffle')
        disp('Shuffle only works with "dont sort" - using descending order')
        sortDir = 'descend';
    end
    mtsorttype = AvAmp;
    [~, sortOrder] = sort(mtsorttype,sortDir);

elseif sortType == 3            % sort by spike amplitude
    if strcmpi(sortDir, 'shuffle')
        disp('Shuffle only works with "dont sort" - using descending order')
        sortDir = 'descend';
    end
    mtsorttype = meanSignal;
    [~, sortOrder] = sort(mtsorttype,sortDir);

elseif sortType == 4            % sort by time of first spike
    if strcmpi(sortDir, 'shuffle')
        disp('Shuffle only works with "dont sort" - using descending order')
        sortDir = 'descend';
    end
    firstSpike(firstSpike == 0) = nan;
    mtsorttype = firstSpike;
    [~, sortOrder] = sort(mtsorttype,sortDir);
end

for i = 1:length(sortOrder)
    inputData(i).sortOrder = i;
    % inputData(sortOrder(i)).filename
end

% sortOrder = [7 9 2 10 3];

processedData = inputData(sortOrder);


processedData(1).riseVector = vertcat(inputData(:).riseTime);
processedData(1).fallVector = vertcat(inputData(:).fallTime);
processedData(1).tauVector = vertcat(inputData(:).tau);
processedData(1).amplitudeVector = vertcat(inputData(:).peakAmplitude);
processedData(1).intervalVector = vertcat(inputData(:).peakIntervals);
processedData(1).meanIntervalVector = vertcat(inputData(:).meanInterval);

if validateRiseFall == 1 && keepValidating == 1
    [file,path] = uiputfile('*.xlsx');
    riseFallSaveName = fullfile(path,file);
    writematrix(processedData(1).riseVector, riseFallSaveName, 'Sheet', 'Rise Time');
    writematrix(processedData(1).fallVector, riseFallSaveName, 'Sheet', 'Fall Time');
    writematrix(processedData(1).tauVector, riseFallSaveName, 'Sheet', 'Decay Constant');
end

if analyzePropagationRate == 1
    processedData(1).propagationVector = vertcat(inputData(:).propagationRate);
    processedData(1).rSquaredVector = vertcat(inputData(:).rSquared);
    processedData(1).validVector = vertcat(inputData(:).validFlags);


    belowThresh = processedData(1).rSquaredVector < settings.minR2;
    aboveThresh = processedData(1).propagationVector;
    aboveThresh(belowThresh) = nan;
    processedData(1).propagationAboveThreshold = aboveThresh;



    if validatePropagationRate == 1
        assignin('base', 'propagationVector', processedData(1).propagationVector);
        assignin('base', 'rSquaredVector', processedData(1).rSquaredVector);
        assignin('base', 'validVector', processedData(1).validVector);

        outputPath = [settings.workingDir '\' processedData(1).genotype '_propagationFit.png'];
        imwrite(fitMontage, outputPath)
    end
end



    disp('Spike Processing Complete')
end

%% Align Traces to Events
function [processedData] = trim2stim(inputData,settings)
for i = 1:length(inputData)
    if isfield(inputData, 'stimTimes')
        if ~isempty(inputData(i).stimTimes)
            firstStim = inputData(i).stimTimes(1);
            prePadding = settings.timePreAlignment*settings.framerate;

            expStart = firstStim-prePadding;
            expEnd = length(inputData(i).bulkSignal);

            if expStart>0
                inputData(i).autoAxialSignal = inputData(i).autoAxialSignal(expStart:expEnd,:);

                inputData(i).bulkSignal = inputData(i).bulkSignal(expStart:expEnd);
                inputData(i).backgroundSignal = inputData(i).backgroundSignal(expStart:expEnd);
                inputData(i).orientation = inputData(i).orientation(expStart:expEnd);
                inputData(i).area = inputData(i).area(expStart:expEnd);
                if isfield(inputData,'velocity')
                    inputData(i).velocity = inputData(i).velocity(expStart:expEnd);
                end



            elseif expStart<1
                nanPadVector = nan(abs(expStart),1);
                nanPadMat = nan(abs(expStart),size(inputData(i).autoAxialSignal,2));

                inputData(i).autoAxialSignal = vertcat(nanPadMat, inputData(i).autoAxialSignal);
                inputData(i).bulkSignal = vertcat(nanPadVector, inputData(i).bulkSignal);
                inputData(i).backgroundSignal = vertcat(nanPadVector,inputData(i).backgroundSignal);
                inputData(i).orientation = vertcat(nanPadVector,inputData(i).orientation);
                inputData(i).area = vertcat(nanPadVector,inputData(i).area);
                if isfield(inputData,'velocity')
                    inputData(i).velocity = vertcat(nanPadVector,inputData(i).velocity);
                end
            end
            inputData(i).stimTimes = inputData(i).stimTimes - expStart;

            if isfield(inputData, 'onFood')
                if ~isempty(inputData(i).onFood)
                    inputData(i).onFood = inputData(i).onFood - expStart;
                end
            end

            if isfield(inputData, 'offFood')
                if ~isempty(inputData(i).offFood)
                    inputData(i).offFood = inputData(i).offFood - expStart;
                end
            end
        end
    end
end
processedData = inputData;
end

%% Process food entrty and exit
function [processedData] = processFood(inputData, settings)
for i = 1:length(inputData)
    if isfield(inputData, 'onFood') && ~isempty(inputData(i).onFood)
        numFrames = find(~isnan(inputData(i).bulkSignal), 1, "last");
        yRange = settings.traceylimit;
        boutData = computeFoodBouts(inputData(i).onFood, inputData(i).offFood, numFrames, yRange);
        inputData(i).onFoodDuration = boutData.onDur/settings.framerate/60;
        inputData(i).offFoodDuration = boutData.offDur/settings.framerate/60;

        inputData(1).onFoodDurationVector = vertcat(inputData(:).onFoodDuration);
        inputData(1).offFoodDurationVector = vertcat(inputData(:).offFoodDuration);
    end
end
processedData = inputData;
end

function [processedData] = alignTraces(inputData,settings)

if settings.trim2stim == 1
    alignmentField = 'stimTimes';
elseif settings.trim2spike == 1
    alignmentField = 'peakLoc';
elseif settings.trim2onFood == 1
    alignmentField = 'onFood';
elseif settings.trim2offFood == 1
    alignmentField = 'offFood';
end


for i = 1:length(inputData)
    
    if isfield(inputData, alignmentField) % make sure the field exists in each dataset
        if settings.trim2stim == 1
            alignmentData = inputData(i).stimTimes;
        elseif settings.trim2spike == 1
            alignmentData = inputData(i).peakLoc;
        elseif settings.trim2onFood == 1
            alignmentData = inputData(i).onFood;
        elseif settings.trim2offFood == 1
            alignmentData = inputData(i).offFood;
            % leavingEvents = inputData(i).offFood;
            % if ~isempty(leavingEvents)
            % 
            %     [priorSpikeIdx] = find(inputData(i).peakLoc < inputData(i).offFood(1), 1, "last");
            % 
            %     if isempty(priorSpikeIdx)
            %         alignmentData = leavingEvents; % trim to leaving event
            %     else
            %         alignmentData = inputData(i).peakLoc(priorSpikeIdx); % trim to last spike before leaving
            %     end
            % end
        end
        
        if ~isempty(alignmentData)
            firstStim = alignmentData(1);
            prePadding = settings.timePreAlignment*settings.framerate;

            expStart = firstStim-prePadding;
            expEnd = length(inputData(i).bulkSignal);

            if expStart>0
                inputData(i).autoAxialSignal = inputData(i).autoAxialSignal(expStart:expEnd,:);

                inputData(i).bulkSignal = inputData(i).bulkSignal(expStart:expEnd);
                inputData(i).backgroundSignal = inputData(i).backgroundSignal(expStart:expEnd);
                inputData(i).orientation = inputData(i).orientation(expStart:expEnd);
                inputData(i).area = inputData(i).area(expStart:expEnd);
                % inputData(i).int1Signal = inputData(i).int1Signal(expStart:expEnd);
                % inputData(i).int9Signal = inputData(i).int9Signal(expStart:expEnd);

                if isfield(inputData,'wormLength')
                    inputData(i).wormLength = inputData(i).wormLength(expStart:expEnd);
                end

                if isfield(inputData,'velocity')
                    inputData(i).velocity = inputData(i).velocity(expStart:expEnd);
                end

                if isfield(inputData,'pumpingRate') && ~isempty(inputData(i).pumpingRate)
                    inputData(i).pumpingRate = inputData(i).pumpingRate(expStart:expEnd);
                end



            elseif expStart<1
                nanPadVector = nan(abs(expStart),1);
                nanPadMat = nan(abs(expStart),size(inputData(i).autoAxialSignal,2));

                inputData(i).autoAxialSignal = vertcat(nanPadMat, inputData(i).autoAxialSignal);
                inputData(i).bulkSignal = vertcat(nanPadVector, inputData(i).bulkSignal);
                inputData(i).backgroundSignal = vertcat(nanPadVector,inputData(i).backgroundSignal);
                inputData(i).orientation = vertcat(nanPadVector,inputData(i).orientation);
                inputData(i).area = vertcat(nanPadVector,inputData(i).area);
                % inputData(i).int1Signal = vertcat(nanPadVector,inputData(i).int1Signal);
                % inputData(i).int9Signal = vertcat(nanPadVector,inputData(i).int9Signal);

                if isfield(inputData,'wormLength')
                    inputData(i).wormLength = vertcat(nanPadVector,inputData(i).wormLength);
                end


                if isfield(inputData,'velocity')
                    inputData(i).velocity = vertcat(nanPadVector,inputData(i).velocity);
                end

                if isfield(inputData,'pumpingRate') && ~isempty(inputData(i).pumpingRate)
                    inputData(i).pumpingRate = vertcat(nanPadVector,inputData(i).pumpingRate);
                end

            end

            if isfield(inputData, 'stimTimes')
                if ~isempty(inputData(i).stimTimes)
                    inputData(i).stimTimes = inputData(i).stimTimes - expStart;
                end
            end

            if isfield(inputData, 'onFood')
                if ~isempty(inputData(i).onFood)
                    inputData(i).onFood = inputData(i).onFood - expStart;
                end
            end

            if isfield(inputData, 'offFood')
                if ~isempty(inputData(i).offFood)
                    inputData(i).offFood = inputData(i).offFood - expStart;
                end
            end
        end
    end
end
processedData = inputData;
end

%% trim experiments
function [processedData] = trimExp(inputData,settings)
lens = nan(length(inputData),1);

for i = 1:length(inputData)
    lens(i)=length(inputData(i).bulkSignal);
end

if settings.trimExperimentLength
    expStart = 1;
    expEnd = min(lens);
elseif settings.analyzePartial
    expStart = settings.partStart;
    expEnd = settings.partEnd;
end

for i = 1:length(inputData)
    inputData(i).autoAxialSignal = inputData(i).autoAxialSignal(expStart:expEnd,:);
    inputData(i).bulkSignal = inputData(i).bulkSignal(expStart:expEnd);
    inputData(i).backgroundSignal = inputData(i).backgroundSignal(expStart:expEnd);
    inputData(i).orientation = inputData(i).orientation(expStart:expEnd);
    inputData(i).area = inputData(i).area(expStart:expEnd);
    % inputData(i).int1Signal = inputData(i).int1Signal(expStart:expEnd);
    % inputData(i).int9Signal = inputData(i).int9Signal(expStart:expEnd);

    if isfield(inputData,'wormLength')
        inputData(i).wormLength = inputData(i).wormLength(expStart:expEnd);
    end

    if isfield(inputData,'velocity')
        inputData(i).velocity = inputData(i).velocity(expStart:expEnd);
    end

    if isfield(inputData,'pumpingRate') && ~isempty(inputData(i).pumpingRate)
        inputData(i).pumpingRate = inputData(i).pumpingRate(expStart:expEnd);
    end

    if isfield(inputData,'stimTimes')
        validtimes = inputData(i).stimTimes<expEnd;
        inputData(i).stimTimes = inputData(i).stimTimes(validtimes);
    end
end
processedData = inputData;
end

%% pad traces with NaN values to automatically equalize exp duration
% will also trim or extend experiments based on the "end" value in analyze partial
% recording (settings.partEnd).
function [processedData] = padTraces(inputData, settings)
lens = nan(length(inputData),1);

for i = 1:length(inputData)
    lens(i)=length(inputData(i).bulkSignal);
end

if settings.partEnd == 0
    maxLen = max(lens);
else
    maxLen = settings.partEnd;
end


for i = 1:length(inputData)
    thisSize = size(inputData(i).autoAxialSignal);
    thisLength = thisSize(1);
    lenDiff = maxLen-thisLength;

    if lenDiff > 0 % Extend traces by padding end with NaN
        traceBuffer = nan(lenDiff,1);
        axialBuffer = nan(lenDiff, thisSize(2));

        inputData(i).autoAxialSignal = vertcat(inputData(i).autoAxialSignal, axialBuffer);
        inputData(i).bulkSignal = vertcat(inputData(i).bulkSignal, traceBuffer);
        inputData(i).backgroundSignal = vertcat(inputData(i).backgroundSignal, traceBuffer);
        inputData(i).orientation = vertcat(inputData(i).orientation, traceBuffer);
        inputData(i).area = vertcat(inputData(i).area, traceBuffer);
        % inputData(i).int1Signal = vertcat(inputData(i).int1Signal, traceBuffer);
        % inputData(i).int9Signal = vertcat(inputData(i).int9Signal, traceBuffer);

        if isfield(inputData,'velocity')
            inputData(i).velocity = vertcat(inputData(i).velocity, traceBuffer);
        end

        if isfield(inputData,'wormLength')
            inputData(i).wormLength = vertcat(inputData(i).wormLength, traceBuffer);
        end

        if isfield(inputData,'pumpingRate') && ~isempty(inputData(i).pumpingRate)
            inputData(i).pumpingRate = vertcat(inputData(i).pumpingRate, traceBuffer);
        end


    elseif lenDiff < 0
        inputData(i).autoAxialSignal = inputData(i).autoAxialSignal(1:maxLen,:);
        inputData(i).bulkSignal = inputData(i).bulkSignal(1:maxLen);
        inputData(i).backgroundSignal = inputData(i).backgroundSignal(1:maxLen);
        inputData(i).orientation = inputData(i).orientation(1:maxLen);
        inputData(i).area = inputData(i).area(1:maxLen);
        % inputData(i).int1Signal = inputData(i).int1Signal(1:maxLen);
        % inputData(i).int9Signal = inputData(i).int9Signal(1:maxLen);

        if isfield(inputData,'velocity')
            inputData(i).velocity = inputData(i).velocity(1:maxLen);
        end

        if isfield(inputData,'wormLength')
            inputData(i).wormLength = inputData(i).wormLength(1:maxLen);
        end

        if isfield(inputData,'pumpingRate') && ~isempty(inputData(i).pumpingRate)
            inputData(i).pumpingRate = inputData(i).pumpingRate(1:maxLen);
        end

        if isfield(inputData,'stimTimes')
            validtimes = inputData(i).stimTimes<maxLen;
            inputData(i).stimTimes = inputData(i).stimTimes(validtimes);
        end


        if isfield(inputData,'onFood')
            validtimes = inputData(i).onFood<maxLen;
            inputData(i).onFood = inputData(i).onFood(validtimes);
        end

        if isfield(inputData,'offFood')
            validtimes = inputData(i).offFood<maxLen;
            inputData(i).offFood = inputData(i).offFood(validtimes);
        end

    end
end
processedData = inputData;
end


function [processedData] = deltaF(inputData, settings)
for i = 1:length(inputData)
    bulkSignal = inputData(i).bulkSignal;
    % %% define F0 as the median of the first 1 second
    % bulkf0 = median(inputData(i).bulkSignal(1:settings.framerate));

    %% define F0 as the median of the lowest 5% of values.
    sortedSignal = sort(bulkSignal,"ascend");
    bulkf0 = mean(sortedSignal(1:round(length(sortedSignal)*0.05)),"omitmissing");
    inputData(i).bulkSignal = sgolayfilt((bulkSignal-bulkf0)./bulkf0,3,25);
end
processedData = inputData;
end




