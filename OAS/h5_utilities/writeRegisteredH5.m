function [] = writeRegisteredH5(path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
d = dir('Z:\OAS\Pha-4GCaMP\250929_pha-4GCaMP_10x_20fps_30gain_noFood2hr\*behavior')

for i = 1:length(d)
    % path = 'C:\Users\Jeremy\Desktop\250723_zfis178_wildtype-multipatch\2025_07_23_21_57_49_flircamera_behavior';
    path = [d(i).folder '\' d(i).name]
    isremote = 1;


    [h5data] = returnRegisteredH5(path, isremote);
    fparts = strsplit(path, '\');




    outputFP = strrep(path, fparts{end}, [fparts{end-1} '_' num2str(i) '.h5']);

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

end