function [behavior_indices, gcamp_indices] = syncDualCameras(beh_folder)
%syncDualCameras takes concatenated timestamps from indexCameraTimestamps()
%and returns a timestamps, file and slice indexes for closest matching
%frames
%
% Usage: Provide a behavior folder as input. Corresponding gcamp recordings
% will be synchronized using the sortest recording as the reference dataset.
% The outputs provide strutures containing the timestamps from the h5
% files "indices.timestamps", an vector indicating which h5 file the frame
% is from "indicies.file_idx", and a vector indicating what slice in that
% file the frame came from "indices.rel_idx". To access a given slice from
% a directory containing multiple h5 files: getslice(folder, file_idx, rel_idx)

if nargin <1
    beh_folder = 'C:\Users\Jeremy\Desktop\250715_zfis178_wildtype-foodEncounter-430_1\2025_07_15_16_49_54_flircamera_behavior';
end

gcamp_folder = strrep(beh_folder, 'behavior', 'gcamp');

[beh_indices] = indexCameraTimestamps(beh_folder);
[gc_indices] = indexCameraTimestamps(gcamp_folder);


%% Plot Camera Jitter
% figure
% plot(beh_t(1:end-1)-beh_t(1), diff(beh_t), gcamp_t(1:end-1)-gcamp_t(1), diff(gcamp_t))
%%

% find the dataset with the fewest values and set this as the reference,
% as we can't register missing frames to a reference. Set the shortest one
% as the reference and the longest as the moving image to register.
if length(beh_indices.timestamps)<length(gc_indices.timestamps)
    register2beh = 1;
    ref_t = beh_indices.timestamps;

    moving_t = gc_indices.timestamps;
    moving_file_idx = gc_indices.h5_file_index;
    moving_rel_idx = gc_indices.relative_index;
else
    register2beh = 0;
    ref_t = gc_indices.timestamps;

    moving_t = beh_indices.timestamps;
    moving_file_idx = beh_indices.h5_file_index;
    moving_rel_idx = beh_indices.relative_index;
end

registered_t = nan(length(ref_t),1);
registered_file_idx = nan(length(ref_t),1);
registered_rel_idx = nan(length(ref_t),1);
frame_mismatch = nan(length(ref_t),1);

% loop through reference dataset and find closeset matching moving dataset
% slice and assign it to the registered dataset
for idx = 1:length(ref_t)
    delta_t = abs(moving_t-ref_t(idx)); % subtract reference from moving time vector
    [frame_diff, match_idx] = min(delta_t); % match_idx is the index of the closest matching frame

    registered_t(idx) = moving_t(match_idx);
    registered_file_idx(idx) = moving_file_idx(match_idx);
    registered_rel_idx(idx) = moving_rel_idx(match_idx);
    frame_mismatch(idx) = frame_diff;
end


% assign the reference and registered timestamps to their appropriate
% camera channel
if register2beh == 1
    behavior_indices = beh_indices;

    gcamp_indices.timestamps = ref_t;
    gcamp_indices.h5_file_index = registered_file_idx;
    gcamp_indices.relative_index = registered_rel_idx;

elseif register2beh == 0
    behavior_indices.timestamps = ref_t;
    behavior_indices.h5_file_index = registered_file_idx;
    behavior_indices.relative_index = registered_rel_idx;

    gcamp_indices = gc_indices;
end

% these arent altered by registration
behavior_indices.h5_path = beh_indices.h5_path;
behavior_indices.img_size = beh_indices.img_size;

gcamp_indices.h5_path = gc_indices.h5_path;
gcamp_indices.img_size = gc_indices.img_size;

