function [] = makeTrackMontage(folder, ax)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% folder  = 'C:\Users\Jeremy\Desktop\2024_1_30_wt_2\2024_01_30_11_30_45_flircamera_behavior';

d = dir([folder '\*videoEvents.mat']);
h5 = dir([folder '\*.h5']);

if ~isempty(d)
    load(fullfile(d(1).folder,d(1).name))
else
    videoEvents = getVideoEvents(folder);
end

x=videoEvents.xSteps;
y=videoEvents.ySteps;

if min(x)<0
    x = (x-min(x))+512;
end

if min(y)<0
    y = (y-min(y))+512;
end

% scatter(x,y)


% img = h5Read(h5file, '/data', [1 1 i], [512,512,1]);
% imshow(img)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example data: x and y coordinates in millimeter units
x_steps = x;
y_steps = y;
% Image resolution and size
% imageResolution = 512;
% pixelSize_mm = 0.002614186;
% imageSize_mm = [1.338463, 1.338463];


% Image resolution and size
imageResolution = 512;
sizeOfOneImage = [512, 512];

% Create a blank canvas for the combined image
imgSz = max(max(y_steps),max(x_steps))+imageResolution;

combinedImageSize = [max(y_steps) + sizeOfOneImage(1) - 1, max(x_steps) + sizeOfOneImage(2) - 1];
combinedImage = zeros(combinedImageSize);

locIdx = length(x_steps);
stepsize = 75;
tic
for j = length(h5):-1:1
    h5file = fullfile(h5(j).folder,h5(j).name);
    info =h5info(h5file, '/data');

    h5Sz = info.Dataspace.Size;
    % Load and position each image onto the combined image
    for i =h5Sz(3):-stepsize:1

        img = fliplr(rot90(imflatfield(h5read(h5file, '/data', [1 1 i], [512,512,1]),20)));
        sizeOfOneImage = [h5Sz(1), h5Sz(2)];  % Image resolution

        % Adjust x and y coordinates based on image size
        if ~isnan(x_steps(locIdx)) && ~isnan(y_steps(locIdx))
            x_position = x_steps(locIdx) - round(sizeOfOneImage(2) / 2); % Adjusted for center alignment
            y_position = y_steps(locIdx) - round(sizeOfOneImage(1) / 2); % Adjusted for center alignment
            % Define the region to update in the combined image
            regionToUpdate = combinedImage(y_position:y_position + sizeOfOneImage(1) - 1, x_position:x_position + sizeOfOneImage(2) - 1);

            % Update the region with non-zero values from the current image
            regionToUpdate(regionToUpdate == 0) = img(regionToUpdate == 0);
%             imshow([regionToUpdate img], [150 256])

            % Update the combined image with the modified region
            combinedImage(y_position:y_position + sizeOfOneImage(1) - 1, x_position:x_position + sizeOfOneImage(2) - 1) = regionToUpdate;
        end
        locIdx = locIdx-stepsize;

        % imshow(combinedImage, [100 256]);
        % % line(x_steps(:),y_steps(:))
        % set(gca, 'YDir', 'normal');
        % drawnow();
    end
end
toc
%% Display the combined image
scatter3(x_steps(:),y_steps(:),1:length(x_steps),1,videoEvents.velocity,Parent=ax);
colormap(ax, turbo)
view(2);
grid off
freezeColors(ax)
hold(ax, 'on');
imh = imshow(combinedImage, [100 250],'Parent',ax);
% colormap bone
hold(ax, 'off')
uistack(imh, 'bottom')

set(gca, 'YDir', 'normal');
ylim([0 combinedImageSize(1)])
xlim([0 combinedImageSize(2)])





end