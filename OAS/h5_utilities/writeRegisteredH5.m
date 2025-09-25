function [] = writeRegisteredH5(path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
path = 'Z:\OAS\unc-31_AID\Defecation\pha-4prom-2_TIR\zfex1268(pha-4prom2-TIR)+4mMAuxPlates\250407_zfis178_UNC-31AID;zfex1268(pha-4-TIR)+4mMAuxPlates_1\2025_04_07_13_57_39_flircamera_behavior';
isremote = 1;


[h5data] = returnRegisteredH5(path, isremote);
fparts = strsplit(path, '\');




outputFP = strrep(path, fparts{end}, [fparts{end-1}, '.h5']);

if isfile(outputFP)
    delete(outputFP)
end

h5create(outputFP, '/gfp', size(h5data.gcamp),'Datatype', 'uint8', 'ChunkSize', size(h5data.gcamp), 'Deflate', 4);
h5create(outputFP, '/bf', size(h5data.brightfield), 'Datatype', 'uint8', 'ChunkSize', size(h5data.brightfield), 'Deflate', 4);
% h5create(outputFP, '/times', size(h5data.time), 'Datatype', 'uint8', 'ChunkSize', size(h5data.time), 'Deflate', 4);
tic
h5write(outputFP, '/gfp',h5data.gcamp)
h5write(outputFP, '/bf', h5data.brightfield);
toc
end