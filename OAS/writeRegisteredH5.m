function [outputArg1,outputArg2] = writeRegisteredH5(behaviorFolder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
behaviorFolder = 'C:\Users\Jeremy\Dropbox\PHB Paper\Supplemental Movies\250522_zfis178_nuc-1-DA-Lactate_5\2025_05_22_14_31_59_flircamera_behavior';

h5data = processH5(behaviorFolder);
fparts = strsplit(behaviorFolder, '\');

outputFP = strrep(behaviorFolder, fparts{end}, [fparts{end-1}, '.h5']);

if isfile(outputFP)
    delete(outputFP)
end

h5create(outputFP, '/gfp', size(h5data.gfp),'Datatype', 'uint8', 'ChunkSize', size(h5data.gfp), 'Deflate', 4);
h5create(outputFP, '/bf', size(h5data.bf), 'Datatype', 'uint8', 'ChunkSize', size(h5data.gfp), 'Deflate', 4);
h5create(outputFP, '/times', size(h5data.time), 'Datatype', 'uint8', 'ChunkSize', size(h5data.time), 'Deflate', 4);
tic
h5write(outputFP, '/gfp',h5data.gfp)
h5write(outputFP, '/bf', h5data.bf);
toc
end