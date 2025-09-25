function [h5data] = returnRegisteredH5(path, isremote)
%returnRegisteredH5 will sync dual channel OAS recordings based on their
%timestamps and return the registered images as a structure
%   Path = path to behavior recording folder
if nargin < 1
    path = 'Z:\OAS\unc-31_AID\Defecation\pha-4prom-2_TIR\zfex1268(pha-4prom2-TIR)+4mMAuxPlates\250407_zfis178_UNC-31AID;zfex1268(pha-4-TIR)+4mMAuxPlates_1\2025_04_07_13_57_39_flircamera_behavior';
    isremote = 1;
end
tic 
[~, nm, ~] = fileparts(path);

behIndPath = ['C:\tmp\' nm 'behInd.mat'];
gcampIndPath = ['C:\tmp\' nm 'gcampInd.mat'];

%% get info from h5 files
if isremote == 1
    if ~isfile(behIndPath) % check if we have any cached data, if not cache it
        disp('Synching Camera feeds...')
        [behavior_indices, gcamp_indices] = syncDualCameras(path,isremote);
        save(behIndPath, 'behavior_indices');
        save(gcampIndPath,'gcamp_indices')
    elseif isfile(behIndPath)           % if we do have cached data, just load it
        disp('loading camera timestamps')
        load(behIndPath)
        load(gcampIndPath)
    end
else
    disp('Synching Camera feeds...')
    [behavior_indices, gcamp_indices] = syncDualCameras(path,isremote);
end


timestamps = behavior_indices.timestamps;
imgWidth = behavior_indices.img_size{1}(1);
imgHeight = behavior_indices.img_size{1}(2);
nFrames = length(timestamps);

brightfield = nan(imgHeight, imgWidth, nFrames);
gcamp = nan(imgHeight, imgWidth, nFrames);

startframe = 1;

for i = startframe:nFrames
    % use synced video indices to check if we need to load the next h5
    % file. Every loop the file index is checked and updates the path
    % of the corresponding h5 file. the relative file index tells which
    % slice to read from that file. When the file index changes, the
    % appropriate h5 file is loaded, otherwise successive slices are
    % read.

    current_beh_file_index = behavior_indices.h5_file_index(i);
    beh_file_path = behavior_indices.h5_path{current_beh_file_index};
    beh_relative_index = behavior_indices.relative_index(i);

    current_gcamp_file_index = gcamp_indices.h5_file_index(i);
    gcamp_file_path = gcamp_indices.h5_path{current_gcamp_file_index};
    gcamp_relative_index = gcamp_indices.relative_index(i);

    %% load the first h5 file
    if i == startframe
        if isremote == 1
            [~, beh_name] = fileparts(beh_file_path);
            [~, gcamp_name] = fileparts(gcamp_file_path);

            beh_local_path = ['C:\tmp\beh' beh_name];
            gcamp_local_path = ['C:\tmp\gcamp' gcamp_name];
            copyfile(beh_file_path, beh_local_path);
            copyfile(gcamp_file_path, gcamp_local_path);

            behavior_h5 = h5read(beh_local_path, '/data');
            gcamp_h5 = h5read(gcamp_local_path, '/data');
            previous_beh_file_index = current_beh_file_index;
            previous_gcamp_file_index = current_gcamp_file_index;

            delete(beh_local_path);
            delete(gcamp_local_path);

        elseif isremote == 0

            behavior_h5 = h5read(beh_file_path, '/data');
            gcamp_h5 = h5read(gcamp_file_path, '/data');
            previous_beh_file_index = current_beh_file_index;
            previous_gcamp_file_index = current_gcamp_file_index;
        end

        %% Load subsequent h5 files when necessary
    elseif i>1
        if current_beh_file_index ~= previous_beh_file_index
            if isremote == 1
                [~, beh_name] = fileparts(beh_file_path);
                beh_local_path = ['C:\tmp\beh' beh_name];
                copyfile(beh_file_path, beh_local_path);
                disp(['loading local copy of: ' beh_file_path])
                behavior_h5 = h5read(beh_local_path, '/data');
                delete(beh_local_path)

            elseif isremote == 0
                disp(['loading ' beh_file_path])
                behavior_h5 = h5read(beh_file_path, '/data');
            end
        end

        if current_gcamp_file_index ~= previous_gcamp_file_index
            if isremote == 1
                [~, gcamp_name] = fileparts(gcamp_file_path);
                gcamp_local_path = ['C:\tmp\beh' gcamp_name];
                copyfile(gcamp_file_path, gcamp_local_path);
                disp(['loading local copy of: ' gcamp_file_path])
                gcamp_h5 = h5read(gcamp_local_path, '/data');
                delete(gcamp_local_path)

            elseif isremote == 0
                disp(['loading ' gcamp_file_path])
                gcamp_h5 = h5read(gcamp_file_path, '/data');
            end
        end


        previous_beh_file_index = current_beh_file_index;
        previous_gcamp_file_index = current_gcamp_file_index;
    end

    

    brightfield(:,:,i) = behavior_h5(:,:, beh_relative_index);
    gcamp(:,:,i) = gcamp_h5(:,:,gcamp_relative_index);

    if mod(i,100) == 0 || i == nFrames
        percentDone = (i / nFrames) * 100;
        fprintf('\rProgress: %3.0f%%', percentDone);
    end
end
h5data.gcamp = gcamp;
h5data.brightfield = brightfield;
h5data.time = timestamps;

toc

end