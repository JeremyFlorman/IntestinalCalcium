function [propagationRate] = getWavePropagationRate(axialPeak, wormLength, settings)


if settings.isOAS == 1
    umPerPixel = 2.6247; % pixel scaling for 4x objective in 2x2 binning on OAS
elseif settings.isOAS == 0
    umPerPixel = 4.9982; % pixel scaling for 5x objective in 4x4 binning on Zeiss inverted scope
end



intestineStart = 30; % where is the boundary between intestine and head (in pixels)
intestineEnd = size(axialPeak,1)-0; % where is the end of the intestine (in pixels)

WormLengthUm = wormLength*umPerPixel;
intestineLengthPx = wormLength*(intestineEnd-intestineStart)/size(axialPeak,1); % Determine # of pixels in the intestine based on how many we deemed part of the head/pharynx
intestineLengthUm = intestineLengthPx*umPerPixel; % convert pixels to micron scaling

numBins = settings.numBins; % how many bins to split the intestine into
binEdges = round(linspace(intestineStart, intestineEnd,numBins+1)); % position of bin edges
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2; % position of bin centers

propMethod =settings.propMethod; % method to define inflection point for propagation rate (1=derivative, 2=half maximum, 3=peak location)
framerate = settings.framerate; % frames per second of recording
validatePropagationRate = settings.validatePropagationRate;
binSizeUm = intestineLengthUm/numBins; % size of bins in microns

fitGlobal = 0; % fit a single linear equation or multiple
numSegments = 2; % if not fitting global equation, how many segments do we want to fit.
fitPad = 0.4; % how much longer to draw fit line

% preallocate variables
binSignal = nan(size(axialPeak,2),numBins);
pkLoc = nan(numBins,1);
pkAmp = nan(numBins, 1);
inflectPt = nan(numBins,1);
fwhm = nan(numBins,1);
changePts = nan(numBins,1);

% x values for plotting binned signal with respect to time
binTime = repmat(linspace(0,size(axialPeak,2)/framerate,size(axialPeak,2))', 1, numBins);

% get average signal in each bin to calculate peak locations and inflection points
for binIdx = 1:numBins
    binSignal(:,binIdx) = mean(axialPeak(binEdges(binIdx):binEdges(binIdx+1),:), 1, 'omitmissing');
    [pk,loc,w] = findpeaks(binSignal(:,binIdx),NPeaks=1, SortStr="descend");

    %% peak location inflection point
    pkLoc(binIdx) = loc; % peak position
    pkAmp(binIdx) = pk; % peak amplitude

    %% derivative inflection point
    % [~, inflectPt(binIdx)] = max(diff(binSignal(:,binIdx))); % find position of the max rate of change in signal
    risePhase = 1:loc;
    dy = smoothdata(gradient(binSignal(risePhase,binIdx)), 'movmean', 5);
    [~, inflectPt(binIdx)] = max(dy);
    %% half maximum peak amplitude inflection point
    fwhm(binIdx) = round(loc-(w/2)); % inflection points defined as points at half-maximum based off half-width from findpeaks
    %% change points
    changePts(binIdx) = findchangepts(binSignal(risePhase,binIdx), "MaxNumChanges",1);
end

waveInit = []; % variable to hold selected inflection points based on what propagation method chosen

if propMethod == 1 % derivative
    waveInit = inflectPt;
elseif propMethod == 2 % full width @ half maximum
    waveInit = fwhm;
elseif propMethod == 3 % peak location
    waveInit = pkLoc;
elseif propMethod == 4
    waveInit = changePts;
end

%% Convert units to microns and seconds
binEdgesUm = binEdges/size(axialPeak,1)*WormLengthUm;
binCentersUm = binCenters/size(axialPeak,1)*WormLengthUm;
initSeconds = waveInit/framerate; % convert frame to seconds

%% Remove Outliers
[~, ~, outlierLoc] = rmoutliers(initSeconds);
% outlierLoc = zeros(size(outlierLoc)); % uncomment to turn off outlier removal

cleanedInit = initSeconds(~outlierLoc);
cleanedCenters = binCentersUm(~outlierLoc);



if ~isempty(cleanedInit)
    %% Calculate Global Propagation Rate
    if fitGlobal == 1
        [coeffs, err]= polyfit(cleanedInit, cleanedCenters, 1); % first value is slope of line fit um/sec
        propagationRate = coeffs(1);

        % plot the polynomial that is used to fit the inflection points
        tFit = linspace(min(cleanedInit)-fitPad, max(cleanedInit)+fitPad, 5);
        yFit = polyval(coeffs, tFit);

        timeValues(1) = {tFit};
        distanceValues(1) = {yFit};
        fitError = err;

    elseif fitGlobal == 0
        %% Calculate Regional Propagation Rate

        regionBoundaries = round(linspace(1, length(cleanedCenters), numSegments+1));
        for i = 1:numSegments
            idxStart = regionBoundaries(i);
            idxEnd = regionBoundaries(i+1)-1;

            xSegment = cleanedInit(idxStart:idxEnd);
            ySegment = cleanedCenters(idxStart:idxEnd);

            % Skip segment if too few points
            if length(xSegment) < 2
                continue
            end

            % Fit line to segment
            [coeffsSegment, err] = polyfit(xSegment, ySegment, 1);
            propagationRate(i) = coeffsSegment(1);
            fitError(i) = err;

            %Optional: Plot fit lines for each segment
            tFit = linspace(min(xSegment)-fitPad, max(xSegment)+fitPad, 5);
            yFit = polyval(coeffsSegment, tFit);

            timeValues(i) = {tFit};
            distanceValues(i) = {yFit};
        end
    end
end


% plotting propagation rates
if validatePropagationRate == 1
    figure(Position=[217.8000 173 447.2000 532.8000], Color=[1 1 1]);
    axFig = tiledlayout(2,1);
    propAx = nexttile();
    axAx = nexttile();
    Col = viridis(numBins);
    colororder(Col);

    %% Trace plotting
    plot(binTime, binSignal,'Parent',propAx)
    hold(propAx, "on")

    for pltIdx =1:numBins
        currentColor = Col(pltIdx,:);
        % plot(pkLoc(pltIdx)/framerate,pkAmp(pltIdx)*1.05, 'v', 'MarkerFaceColor', currentColor, 'MarkerEdgeColor', currentColor, 'Parent', propAx);
        if ~isempty(initSeconds(pltIdx))
            initPt = initSeconds(pltIdx);
            plot(initPt, binSignal(waveInit(pltIdx),pltIdx)*1.05, 'v', 'MarkerFaceColor', currentColor, 'MarkerEdgeColor', currentColor,'Parent',propAx)
        end
    end

    hold(propAx, "off")
    xlim(propAx, [0 binTime(end)])


    %% Kymograph plotting
    imagesc(binTime(:,1), linspace(0,WormLengthUm,size(axialPeak,1)), axialPeak,'Parent', axAx)
    colormap(axAx,'turbo');

    for pltIdx = 1:numBins
        rectX = 0;
        rectY = binEdgesUm(pltIdx); % this should be pltIdx+1 for the lower left corner of rectangle, but Y axis is reversed!
        rectW = size(binSignal,1);
        rectH = binSizeUm;

        currentColor = Col(pltIdx,:);

        % plot rectangle around each bin
        % rectangle('Position',[rectX, rectY,rectW, rectH],'linestyle', '-','EdgeColor', currentColor, 'linewidth', 1, 'Parent', axAx)

        % plot a line at each inflection point
        if outlierLoc(pltIdx) == 0
            line([initSeconds(pltIdx) initSeconds(pltIdx)], [rectY rectY+rectH], 'Color', [1 1 1], 'linewidth' ,1.5, 'Parent', axAx)
        elseif outlierLoc(pltIdx) == 1
            line([initSeconds(pltIdx) initSeconds(pltIdx)], [rectY rectY+rectH], 'Color', [1 0 0], 'linewidth' ,1.5, 'Parent', axAx)
        end
    end



    %% Plot Lines Fit to Each Segment

    for i=1:length(timeValues)
        line(timeValues{i}, distanceValues{i}, 'Color', 'r', 'Linestyle', ':', 'LineWidth', 1.5);
    end

    %% Axis Labels and Formatting
    ylabel(propAx, 'GCaMP Signal (a.u)')
    xlabel(propAx, 'Time (sec)')
    box(propAx, 'off')

    propAx.YAxis.Exponent = 4;

    ylabel(axAx, 'Distance (\it{\mum})')
    xlabel(axAx, 'Time (sec)')
    box(axAx, 'off')

    cb = colorbar(axAx);
    cb.Ruler.Exponent = 4;


    title(['Slope: ' num2str(round(propagationRate,1)) ' \mum/sec'] ,  ['R^2 ' num2str([fitError(:).rsquared])], 'Parent',axFig)
    % txt = input("Look ok? if not press letter key before hitting enter. Type 'exit' to quit","s");
    %
    % if ~isempty(txt)
    %     propagationRate = NaN;
    % end

end
end
