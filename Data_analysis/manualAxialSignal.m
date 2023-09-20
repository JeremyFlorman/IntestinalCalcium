function [] = manualAxialSignal(folder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

folder = uigetdir('X:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\TA_resistance\lgc-55-int-resc-1-30mM-TA\');

textfilepath = dir([folder '\*reslice*']);
matpath = dir([folder '\*wormdata.mat']);

load(fullfile(matpath.folder, matpath.name));

data = readmatrix(fullfile(textfilepath.folder,textfilepath.name));
w = size(wormdata.autoAxialSignal,2);

for i = 1:length(data)
    wormdata.autoAxialSignal(i,1:w) = resample(data(i,:), w, size(data,2),5,20);
end

save(fullfile(matpath.folder, matpath.name),"wormdata");

imagesc(wormdata.autoAxialSignal)