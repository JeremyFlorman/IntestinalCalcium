function [indices] = indexCameraTimestamps(folder)
%indexCameraTimestamps concatenates timestamps in chunked video recordings
%and returns the times, the cooresponding h5 file and the relative index
%within the file


if nargin <1
    folder = 'C:\Users\Jeremy\Desktop\250715_zfis178_wildtype-foodEncounter-430_1\2025_07_15_16_49_54_flircamera_behavior';
end

h5_files = dir([folder '\*.h5']);

timestamps = [];
h5_file_index = [];
relative_index = [];
h5_path = cell(3,1);
img_size = cell(3,1);

for idx = 1:length(h5_files)
    disp(['Getting timestamps from file ' num2str(idx) ' of ' num2str(length(h5_files))])
    this_file = fullfile(h5_files(idx).folder, h5_files(idx).name);
    local_filepath = ['C:\tmp\' h5_files(idx).name];
    copyfile(this_file, local_filepath);

    temp_times = h5read(this_file, '/times');
    temp_relative_index = 1:size(temp_times,1);
    timestamps = vertcat(timestamps, temp_times);
    h5_file_index = vertcat(h5_file_index, repmat(idx, size(temp_times,1), 1));
    relative_index = vertcat(relative_index, temp_relative_index');
    h5_path(idx) = {this_file};

    img_info = h5info(this_file, '/data');
    img_size(idx) = {img_info.ChunkSize};

    delete(local_filepath);
end


indices.timestamps = timestamps;
indices.h5_file_index = h5_file_index;
indices.relative_index = relative_index;
indices.h5_path = h5_path;
indices.img_size = img_size;