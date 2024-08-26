function [] = drawSkeletonROI()
% Pulls up a frame and allows the user to draw an ROI along the worm spine
% for use in extracting axial signal. This can be helpful for situations
% when worms coil and remain stationary. 
%   Detailed explanation goes here


path = 'Z:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\TA_resistance\lgc-55-30mM-TA\220823_zfis178_lgc-55-30mM-TA_1';
frame2Read = 1; 
cropPixels = 3;

t = dir([path '\*.tif']);

img=imread(fullfile(t(1).folder,t(1).name), frame2Read);

img = img(cropPixels:end-cropPixels,cropPixels:end-cropPixels);

figure();
imshow(img, [3000 25000])

roi = drawpolyline();

outname = strrep(fullfile(t(1).folder,t(1).name), 'MMStack_Default.ome.tif', ['roi_' num2str(ceil(frame2Read/2)) '.mat']);
save(outname, 'roi')
end