function [outputArg1,outputArg2] = alignOASVideo(bPath, gPath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bPath = "C:\Users\Jeremy\Desktop\ZFIS178\2023_07_21_15_26_29_flircamera_behavior\000000.h5";
gPath = "C:\Users\Jeremy\Desktop\ZFIS178\2023_07_21_15_26_29_flircamera_gcamp\000000.h5";

bData = h5read(bPath, '/data');
gData = h5read(gPath, '/data');

btime = h5read(bPath, '/times');
gtime = h5read(gPath, '/times');

[imgWidth, imageHeight, nFrames] = size(gData);

plot(btime,gtime)


end