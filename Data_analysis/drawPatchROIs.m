folder  = 'C:\Users\Jeremy\Desktop\260325_zfis178_wildtype-7patch-r12mm-05ul-15fps_1\2026_03_25_11_54_09_flircamera_behavior';
tic
d  = dir(fullfile(folder, '*videoEvents.mat'));
h5 = dir(fullfile(folder, '*.h5'));

if ~isempty(d)
    load(fullfile(d(1).folder, d(1).name));   % loads videoEvents
else
    videoEvents = getVideoEvents(folder);
end

%% 1) Use mm coordinates directly
x_mm = videoEvents.xLoc;   % in mm
y_mm = videoEvents.yLoc;   % in mm

% Make sure these are column vectors
x_mm = x_mm(:);
y_mm = y_mm(:);

% Shift so everything is positive (preserves relative geometry)
x_mm = x_mm - min(x_mm);
y_mm = y_mm - min(y_mm);

%% 2) Convert mm → pixels
pixelSize_mm = 0.002614186;
pxPerMm      = 1 / pixelSize_mm;        % ≈ 382.53 px / mm

x_px_center = x_mm * pxPerMm;
y_px_center = y_mm * pxPerMm;

% Round once to pixel centers
x_px_center = round(x_px_center);
y_px_center = round(y_px_center);

imageResolution = 512;                   % 512x512
halfSize = imageResolution / 2;

nEvents = numel(x_px_center);

%% 3) Determine canvas size in pixels
% Rough bounds for centers
minXc = min(x_px_center);
minYc = min(y_px_center);
maxXc = max(x_px_center);
maxYc = max(y_px_center);

% Account for half the image on each side
minX = floor(minXc - halfSize);
minY = floor(minYc - halfSize);
maxX = ceil(maxXc + halfSize);
maxY = ceil(maxYc + halfSize);

% Shift offsets so indices are 1-based
offsetX = 1 - minX;
offsetY = 1 - minY;

canvasWidth  = maxX - minX + 1;
canvasHeight = maxY - minY + 1;

combinedImage = inf(canvasHeight, canvasWidth);  % start "very bright" (for min projection)

%% 4) Walk through H5s and map frames to positions by global frame index
% (Assumes videoEvents has one row per frame, in time order)
% Sort h5 files by name just to be safe
[~, idx] = sort({h5.name});
h5 = h5(idx);

frameOffset = 0;
stepSize    = 75;   %Subsampling

for j = 1:numel(h5)
    h5file = fullfile(h5(j).folder, h5(j).name);
    info   = h5info(h5file, '/data');
    h5Sz   = info.Dataspace.Size;   % [height width nFrames]
    nFrames = h5Sz(3);

    for f = 1:stepSize:nFrames
        k = frameOffset + f;       % global frame index

        if k > nEvents
            break;                 % no more coordinates
        end

        if ~isfinite(x_px_center(k)) || ~isfinite(y_px_center(k))
            continue;
        end

        % 5) Read this frame
        rawFrame = h5read(h5file, '/data', [1 1 f], [h5Sz(1) h5Sz(2) 1]);
        rawFrame = double(rawFrame);           % to match combinedImage
        img      = rot90(rawFrame);  % your flat-field correction

        % 6) Compute top-left placement in canvas coords
        xc = x_px_center(k) + offsetX;
        yc = y_px_center(k) + offsetY;

        x1 = round(xc - halfSize + 1);
        y1 = round(yc - halfSize + 1);
        x2 = x1 + imageResolution - 1;
        y2 = y1 + imageResolution - 1;

        % sanity check to keep in bounds
        if x1 < 1 || y1 < 1 || x2 > canvasWidth || y2 > canvasHeight
            % If this ever happens a lot, we can pad or adjust offsets,
            % but for now just skip out-of-bounds frames.
            continue;
        end

        % 7) Update background using per-pixel minimum
        region = combinedImage(y1:y2, x1:x2);   % existing region
        combinedImage(y1:y2, x1:x2) = min(region, img);  % darker wins
    end

    frameOffset = frameOffset + nFrames;
end

% Replace infinities (never written pixels) with something reasonable
% e.g. median of all finite pixels:
finiteVals = combinedImage(isfinite(combinedImage));
bgFill     = median(finiteVals, 'omitnan');
combinedImage(~isfinite(combinedImage)) = bgFill;


%% 5) Visualize with trajectory overlay (using same pixel coords)
[pth, ~, ~] = fileparts(folder);
matdir = dir([pth '\*wormdata.mat' ]);

% save combined image at full resolution w/o overlays
% img_disp = mat2gray(combinedImage, [prctile(finiteVals,1)  prctile(finiteVals,99)]);
% imwrite(img_disp, strrep(fullfile(matdir.folder, matdir.name),'wormdata.mat','combinedImage.png'));
if isscalar(matdir)
    wormDataPath = fullfile(matdir(1).folder, matdir(1).name);
    load(wormDataPath)

else
    [dataFile, dataPath]= uigetfile([pth '\*wormdata.mat' ]);
    wormDataPath = fullfile(dataPath, dataFile);
    load(wormDataPath)
end

peakLoc = wormdata.peakLoc;

bulkSignal = smoothdata(wormdata.bulkSignal-wormdata.backgroundSignal, 'movmean', 60);

[int1Signal, int9Signal] = extractCellSignal(wormdata); % Get cell specific signal
int1Signal = int1Signal - wormdata.backgroundSignal;    % subtract background
int9Signal = int9Signal - wormdata.backgroundSignal;
inc   = 1:75:nEvents; % subsample for faster plotting

%% Correct head/tail swaps
if isfield(wormdata, 'headLoc') && isfield(wormdata, 'tailLoc')
    headXY = wormdata.headLoc;
    tailXY = wormdata.tailLoc;

    nFrames = size(headXY,1);

    cleanedHead = headXY;
    cleanedTail = tailXY;

    swapped = false(nFrames,1);

    for i = 2:nFrames
        prevHead = cleanedHead(i-1,:);
        prevTail = cleanedTail(i-1,:);

        currHead = headXY(i,:);
        currTail = tailXY(i,:);

        % cost if current labels are correct
        cost_noSwap = norm(currHead - prevHead) + norm(currTail - prevTail);

        % cost if current labels should be swapped
        cost_swap = norm(currTail - prevHead) + norm(currHead - prevTail);

        if cost_swap < cost_noSwap
            cleanedHead(i,:) = currTail;
            cleanedTail(i,:) = currHead;
            swapped(i) = true;
        else
            cleanedHead(i,:) = currHead;
            cleanedTail(i,:) = currTail;
        end
    end


    %% convert Head and Tail coordinates to common frame of reference
    xCenter = x_px_center;
    yCenter = y_px_center;
    frameCenter = [xCenter yCenter];

    imgCenter = [(imageResolution+1)/2, (imageResolution+1)/2];   % matches stitch placement convention

    % Correct for 90 degree rotation
    headLocal = [cleanedHead(:,2), imageResolution + 1 - cleanedHead(:,1)];
    tailLocal = [cleanedTail(:,2), imageResolution + 1 - cleanedTail(:,1)];

    % fill missing values
    headGlobal = fillmissing(frameCenter + (headLocal - imgCenter) + offsetX, 'next');
    midGlobal = fillmissing([xCenter+offsetX yCenter+offsetY], 'linear');
    tailGlobal = fillmissing(frameCenter + (tailLocal - imgCenter) + offsetY, 'next');

    % Subsample for plotting
    xSubset = headGlobal(inc,1);
    ySubset = headGlobal(inc,2);
else

    midGlobal = fillmissing([xCenter+offsetX yCenter+offsetY], 'linear');

    % Subsample for plotting
    xSubset = x_px_center(inc) + offsetX;
    ySubset = y_px_center(inc) + offsetY;
end
%%
colorSignal = int9Signal; % signal used for color of scatter plot markers
sizeSignal = repmat(8, size(colorSignal)); % signal used for size of scatter plot markers
% sizeSignal(sizeSignal<5) = 1;


figure;
ax = gca;


s = scatter(ax,xSubset,ySubset,sizeSignal(inc),colorSignal(inc), 'o', 'filled');
colormap(ax, turbo);
ax.CLim = [0 45];

view(2);
grid off
freezeColors(ax)
hold(ax, 'on');

% scatter(x_px_center(peakLoc)+offsetX, y_px_center(peakLoc)+offsetY+150, 20, 'yv','filled', 'MarkerEdgeColor', 'k')


% cb = colorbar;
% cb.Label.String = 'GCaMP Signal (a.u.)';
% cb.Label.Rotation = -90;
imh = imshow(combinedImage, [prctile(finiteVals,1) prctile(finiteVals,99)]);
uistack(imh, 'bottom');


set(ax, 'YDir', 'normal');
ylim(ax, [1 canvasHeight]);
xlim(ax, [1 canvasWidth]);


% Convert axis labels from pixels to mm
mmTickStep = 2;
xlim_mm = ax.XLim / pxPerMm;
ylim_mm = ax.YLim / pxPerMm;

% Choose round-number mm ticks within the visible range
xticks_mm = ceil(xlim_mm(1)/mmTickStep)*mmTickStep : ...
    mmTickStep : ...
    floor(xlim_mm(2)/mmTickStep)*mmTickStep;

yticks_mm = ceil(ylim_mm(1)/mmTickStep)*mmTickStep : ...
    mmTickStep : ...
    floor(ylim_mm(2)/mmTickStep)*mmTickStep;

% Convert back to pixels for positioning
xticks_px = xticks_mm * pxPerMm;
yticks_px = yticks_mm * pxPerMm;

set(gca, ...
    'XTick', xticks_px, ...
    'YTick', yticks_px, ...
    'XTickLabel', arrayfun(@num2str, xticks_mm, 'UniformOutput', false), ...
    'YTickLabel', arrayfun(@num2str, yticks_mm, 'UniformOutput', false));

xlabel('');
ylabel('');

ax.XAxis.Visible = 0;
ax.YAxis.Visible = 0;

hold(ax, 'off');
toc
%% Define food patches

% ROIs = struct();
nPatches = input("How Many Food Patches?");
if ~isempty(nPatches) && isnumeric(nPatches)
    for i = 1:nPatches
        disp(['Draw Circle ' num2str(i) ' of ' num2str(nPatches) ...
            ' - Double Click When Done'])
        c = drawcircle('Color',[0.6350 0.0780 0.1840], FaceAlpha=0, LineWidth=1);
        wait(c)
        ROIs(i).Center = c.Center;
        ROIs(i).Radius = c.Radius;
        ROIs(i).Vertices = c.Vertices;
    end
end

%% Set up query points for ROI
if isfield(wormdata, 'headLoc') && isfield(wormdata, 'tailLoc')
    xq= headGlobal(:,1); % define x querry points
    yq = headGlobal(:,2); % define y querry points

ROIs(1).headXY = headGlobal;
ROIs(1).midXY = midGlobal;
ROIs(1).tailXY = tailGlobal;

nPoints = numel(xq);
headIn = nan(nPoints, nPatches);
midIn = nan(nPoints, nPatches);
tailIn = nan(nPoints, nPatches);
else
    xq= midGlobal(:,1); % define x querry points
    yq = midGlobal(:,2); % define y querry points
    ROIs(1).midXY = midGlobal;
    
    nPoints = numel(xq);
    inPoints = nan(nPoints, nPatches);
end



% Compute points inside patches
for i =1:numel(ROIs)

    c = drawcircle('Center', ROIs(i).Center,'Radius', ...
        ROIs(i).Radius, 'Color',[0.9059    0.1608    0.5412], ...
        FaceAlpha=0, LineWidth=1.5, MarkerSize=1); % dilate circle by 1mm

    if isfield(wormdata, 'headLoc') && isfield(wormdata, 'tailLoc')
        headIn(1:nPoints, i) = inROI(c, headGlobal(:,1), headGlobal(:,2));
        midIn(1:nPoints, i) = inROI(c, midGlobal(:,1), midGlobal(:,2));
        tailIn(1:nPoints, i) = inROI(c, tailGlobal(:,1), tailGlobal(:,2));
    else 
        inPoints(1:nPoints, i) = inROI(c, midGlobal(:,1), midGlobal(:,2));
    end
end


if isfield(wormdata, 'headLoc') && isfield(wormdata, 'tailLoc')
    allIn = any(headIn,2) | any(midIn,2); % Animal is on food if its head or midbody are on food,
                                           % this will exclude headpokes outside of the foodpatch from being considered leaving events
else 
    allIn = any(inPoints, 2);
end



wormdata.onFoodVector = allIn;
wormdata.patchROIs = ROIs;

boutData = computeFoodBouts(allIn, 15, 0);
wormdata.boutData = boutData;

save(wormDataPath,"wormdata")
clear('ROIs')


plot_SummaryTraces(wormDataPath)

pc = computePhaseChange(wormdata,15, folder);



% figure()
% plot(xq(allIn), yq(allIn), 'yo')
% hold on
% plot(xq(~allIn), yq(~allIn), 'rx')
%%
%
%
%
% [int1Clipped, ~] = clipSpikes(wormdata.peakLoc, 750, int1Signal);
% [int9Clipped, ~] = clipSpikes(wormdata.peakLoc, 750, int9Signal);
%
% toClip = 1;
% if toClip == 1
%     int1Plotting = int1Clipped;
%     int9Plotting = int9Clipped;
% else
%     int1Plotting = int1Signal;
%     int9Plotting = int9Signal;
% end
%
%
%
% int1onFood = int1Plotting;
% int1onFood(~allIn) = nan;
%
% int1offFood = int1Plotting;
% int1offFood(allIn) = nan;
%
% int9onFood = int9Plotting;
% int9onFood(~allIn) = nan;
%
% int9offFood = int9Plotting;
% int9offFood(allIn) = nan;
%
% time = linspace(0, length(bulkSignal)/30/60, length(bulkSignal));
%
% [onFood, offFood] = transitionPts(allIn);
%
%
% % annotate food patches
% xpatch = nan(4, length(onFood));
% for i=1:length(onFood)
%     s = onFood(i,1);
%     e = onFood(i,2);
%     xpatch(1:4, i) = [s s e e];
% end
% xpatch = xpatch/30/60;
% yhi = max([int1Plotting int9Plotting], [],'all', 'omitmissing')*1.1;
% ypatch = repmat([0; yhi; yhi; 0],1 ,length(onFood));
%
% if toClip == 0
%     histYLims = [0 0.25];
%     histXLims = [0 80];
%     plotLims = [0 yhi];
% elseif toClip == 1
%     histYLims = [0 0.25];
%     histXLims = [0 80];
%     plotLims = [0 yhi];
% end
% histBins = 0:2:100;
%
% figure(Color=[1 1 1], ...
%     Position=[2397 -3.8000 1.1104e+03 392.8000]);
%
% tl = tiledlayout(2,3, 'TileSpacing','tight', 'Padding','compact');
%
% % Int1 Signal
% nexttile([1 2])
% plot(time, int1Plotting, 'k')
%
%
% p = patch(xpatch, ypatch, [0.93 0.69 0.13], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
% uistack(p, 'bottom');
% ylim(plotLims)
% ylabel('Int1 Signal');
% xlabel('Time (min)')
% box off
%
%
% nexttile([1 1])
% h1 = histogram(int1onFood,'BinEdges', histBins, 'Normalization','probability','EdgeAlpha',0.25);
% hold on
% h2 = histogram(int1offFood,'BinEdges', histBins,'Normalization','probability','EdgeAlpha',0.25);
% hold off
% legend();
%
% ylim(histYLims)
% xlim(histXLims)
% ylabel('Probability');
% xlabel('Int1 Signal')
% box off
%
%
% nexttile([1 2])
% plot(time, int9Plotting, 'k')
% p = patch(xpatch, ypatch, [0.93 0.69 0.13], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
% uistack(p, 'bottom');
% ylim(plotLims)
% ylabel('Int9 Signal');
% xlabel('Time (min)')
% box off
%
% nexttile([1 1])
%
% h3 =histogram(int9onFood,'BinEdges', histBins,  'Normalization','probability', 'EdgeAlpha',0.25);
% hold on
% h4 =histogram(int9offFood,'BinEdges', histBins,  'Normalization','probability','EdgeAlpha',0.25);
% hold off
% legend();
% ylim(histYLims)
% xlim(histXLims)
% ylabel('Probability');
% xlabel('Int9 Signal')
% box off
% %%
%
%
% % [onFood offFood]
%
%
% % plot(time, bulkSignal, time, clippedSignal)
%
% function [clippedSignal, spikeFlags] = clipSpikes(spikes, windowFrames, signal)
%
% spikeFlags = true(length(signal), 1);
%
% for i =1:length(spikes)
%     preSpikeIdx = spikes(i)-floor(windowFrames*0.33);
%     postSpikeIdx = spikes(i)+floor(windowFrames*0.66);
%
%     spikeFlags(preSpikeIdx:postSpikeIdx, 1) = 0;
% end
%
% clippedSignal = signal;
% clippedSignal(~spikeFlags) = nan;
% end
%
%
%
%
% function [onBouts, offBouts] = transitionPts(foodVector)
% onFood = [];
% offFood = [];
%
% %% find whether the worm starts on or off food
% switch foodVector(1)
%     case 1
%         onFood(1) = 1;
%     case 0
%         offFood(1) = 1;
% end
%
% %% find transition points in the foodVector (1 marks on food, 0 off food)
% dif = diff(foodVector); % determine transitions by subracting the previous value, 0 = no change
% off = find(dif==-1); % -1 is a leaving event
% on = find(dif==1); % +1 is an entry event
%
% %% concatenate the starting state
% onFood = vertcat(onFood, on);
% offFood = vertcat(offFood, off);
%
%
% %% compute 'on food' bouts (beginning and end frames)
% onBouts = nan(length(onFood),2);
%
% for i =1:length(onFood)
%     thisOn = onFood(i);
%     onBouts(i,1) = thisOn;
%
%     % find the subsequent leaving event
%     thisOff = find(offFood>thisOn, 1, "first");
%
%     if ~isempty(thisOff) % make sure its not the last timepoint
%         onBouts(i,2) = offFood(thisOff);
%     else
%         onBouts(i,2) = length(foodVector);
%     end
% end
%
%
% %% compute 'off food' bouts (beginning and end frames)
% offBouts = nan(length(offFood),2);
%
% for i =1:length(offFood)
%     thisOff = offFood(i);
%     offBouts(i,1) = thisOff;
%
%     % find the subsequent entry event
%     thisOn = find(onFood>thisOff, 1, "first");
%     if ~isempty(thisOn) % make sure its not the last timepoint
%         offBouts(i,2) = onFood(thisOn);
%     else
%         offBouts(i,2) = length(foodVector);
%     end
% end
%
% end
%
%
%
