folder  = 'Z:\OAS\foodEncounter\wildtype-noFood-0minCycleStart\260209_zfis178_wildtype-noFood-0min_1\2026_02_09_14_09_10_flircamera_behavior';

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
halfSize       = imageResolution / 2;

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
stepSize    = 75;   % your subsampling

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


%% 8) Visualize with trajectory overlay (using same pixel coords)
[pth, ~, ~] = fileparts(folder);
matdir = dir([pth '\*wormdata.mat' ]);

% save combined image at full resolution w/o overlays
% img_disp = mat2gray(combinedImage, [prctile(finiteVals,1)  prctile(finiteVals,99)]);
% imwrite(img_disp, strrep(fullfile(matdir.folder, matdir.name),'wormdata.mat','combinedImage.png'));

load(fullfile(matdir(1).folder, matdir(1).name))
peakLoc = wormdata.peakLoc;

bulkSignal = smoothdata(wormdata.bulkSignal-wormdata.backgroundSignal, 'movmean', 60);

[int1Signal, int9Signal] = extractCellSignal(wormdata); % Get cell specific signal
int1Signal = int1Signal - wormdata.backgroundSignal;    % subtract background
int9Signal = int9Signal - wormdata.backgroundSignal;


colorSignal = int1Signal; % signal used for color of scatter plot markers
sizeSignal = abs(int1Signal-5); % signal used for size of scatter plot markers
inc   = 1:5:nEvents; % subsample for faster plotting

figure;
ax = gca;
s = scatter(ax, x_px_center(inc) + offsetX, y_px_center(inc) + offsetY,sizeSignal(inc),colorSignal(inc), 'o', 'filled');
colormap(ax, turbo);
ax.CLim = [0 60];

view(2);
grid off
freezeColors(ax)
hold(ax, 'on');

scatter(x_px_center(peakLoc)+offsetX, y_px_center(peakLoc)+offsetY+150, 20, 'yv','filled', 'MarkerEdgeColor', 'k')


% cb = colorbar;
% cb.Label.String = 'Time (min)';
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

% xlabel('mm');
% ylabel('mm');

% ax.XAxis.Visible = 1;
% ax.YAxis.Visible = 1;

hold(ax, 'off');
