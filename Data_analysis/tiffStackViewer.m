function tiffStackViewer(startpath)
% Create the figure and UI components

if isempty(startpath)
    % Load TIFF stack
    [file, path] = uigetfile([startpath '\*.tif'], 'Select a multi-slice TIFF file');
    if isequal(file, 0)
        disp('User selected Cancel');
        return;
    else
        fullFilePath = fullfile(path, file);
    end
else
    d = dir([startpath '\*.tif']);
    path = d(1).folder;
    file = d(1).name;
    fullFilePath = fullfile(path, file);
end




% copy file locally
tempfolder = 'C:\tmp';
localpath = fullfile(tempfolder,file);
if ~exist(localpath, "file")
    tic
    disp(['Copying file to: ' tempfolder ' for local processing'])
    copyfile(fullFilePath, tempfolder)
    disp(['File copied in ' num2str(toc)])
end

% Initialize variables

tiff_info = imfinfo(localpath);
num_slices = numel(tiff_info);
cropPixels = 3;
current_frame = 1;


% load existing ROIs if they exist
roiname = strrep(fullFilePath, 'MMStack_Default.ome.tif', 'ROIs.mat');
if exist(roiname,"file")
    load(roiname);
    disp('ROIs loaded successfully!');
else
    disp('No previous ROIs found, starting from scratch...')
    roiStruct = struct('roi', {}, 'frame', {});
end


hFig = figure('Name', 'TIFF Stack Viewer', 'NumberTitle', 'off', ...
    'Position', [100, 100, 600, 400]);

% Store variables in the figure's UserData
set(hFig, 'UserData', struct('roi', [], 'roiStruct', roiStruct, 'current_frame', current_frame));



% Create image display axis
hAx = axes('Parent', hFig, 'Units', 'pixels', 'Position', [50, 100, 500, 250]);
img = imread(localpath, current_frame);
img = img(cropPixels:end-cropPixels, cropPixels:end-cropPixels);
hImg = imshow(img, [], 'Parent', hAx);

% Create slider for frame selection
hSlider = uicontrol('Style', 'slider', 'Min', 1, 'Max', num_slices, ...
    'Value', 1, 'SliderStep', [1/(ceil(num_slices/2)-1), 1/(ceil(num_slices/2)-1)], ...
    'Position', [100, 50, 400, 20], ...
    'Callback', @slider_callback);

% Create frame number entry box
hFrameEdit = uicontrol('Style', 'edit', 'String', '1', ...
    'Position', [510, 51.2, 80, 20], ...
    'Callback', @edit_callback);

% Create "Draw Spine" button
hDraw = uicontrol('Style', 'pushbutton', 'String', 'Draw Spine', ...
    'Position', [150, 28, 80, 20], ...
    'Callback', @draw_callback);

% Create "Add to structure" button
hAdd = uicontrol('Style', 'pushbutton', 'String', 'Add to structure', ...
    'Position', [250, 28, 80, 20], ...
    'Callback', @add_callback);

% Create "AddStopp to structure" button
hAddStop = uicontrol('Style', 'pushbutton', 'String', 'Add a Stop', ...
    'Position', [250, 5, 80, 20], ...
    'Callback', @addStop_callback);


% Create "Restore ROI" button
hRestore = uicontrol('Style', 'pushbutton', 'String', 'Restore ROI', ...
    'Position', [350, 28, 80, 20], ...
    'Callback', @restore_callback);

% Create "Save ROIs" button
hSaveButton = uicontrol('Style', 'pushbutton', 'String', 'Save ROIs', ...
    'Position', [450, 28, 80, 20], ...
    'Callback', @save_callback);


% Slider callback function
    function slider_callback(hObject, ~)
        userData = get(hFig, 'UserData');
        userData.current_frame = round(get(hObject, 'Value'));
        set(hFig, 'UserData', userData);
        update_image();
    end

% Edit box callback function
    function edit_callback(hObject, ~)
        frame_num = str2double(get(hObject, 'String'));
        userData = get(hFig, 'UserData');
        if isnan(frame_num) || frame_num < 1 || frame_num > num_slices
            set(hObject, 'String', num2str(userData.current_frame));
            return;
        end
        userData.current_frame = round(frame_num);
        set(hSlider, 'Value', userData.current_frame);
        set(hFig, 'UserData', userData);
        update_image();
    end

% Draw callback function
    function draw_callback(~, ~)
        userData = get(hFig, 'UserData');
        if ~isempty(userData.roi)
            delete(userData.roi); % Delete the previous ROI if exists
        end
        userData.roi = drawpolyline(); % Draw a new ROI
        set(hFig, 'UserData', userData);
    end

% Add to structure callback function
    function add_callback(~, ~)
        userData = get(hFig, 'UserData');
        if isempty(userData.roi)
            warndlg('No ROI to add. Please draw or restore an ROI first.');
            return;
        end

        % Add the new ROI and frame information to the structure
        newRoiData = struct('roi', userData.roi.Position, 'frame', userData.current_frame);

        % Append the new ROI data to the roiStruct
        userData.roiStruct(end+1) = newRoiData;

        % Sort the roiStruct by frame to maintain order
        [~, sortIdx] = sort([userData.roiStruct.frame]);
        userData.roiStruct = userData.roiStruct(sortIdx);

        % Clear the ROI from the image
        delete(userData.roi);
        userData.roi = []; % Clear the ROI from UserData
        set(hFig, 'UserData', userData);
    end

% Add to structure callback function
    function addStop_callback(~, ~)
        userData = get(hFig, 'UserData');
        % if isempty(userData.roi)
        %     warndlg('No ROI to add. Please draw or restore an ROI first.');
        %     return;
        % end

        % Add the new ROI and frame information to the structure
        newRoiData = struct('roi', [], 'frame', userData.current_frame);

        % Append the new ROI data to the roiStruct
        userData.roiStruct(end+1) = newRoiData;

        % Sort the roiStruct by frame to maintain order
        [~, sortIdx] = sort([userData.roiStruct.frame]);
        userData.roiStruct = userData.roiStruct(sortIdx);

        set(hFig, 'UserData', userData);
    end


% Restore ROI callback function
    function restore_callback(~, ~)
        userData = get(hFig, 'UserData');
        roiStruct = userData.roiStruct;
        if isempty(roiStruct)
            warndlg('No ROIs have been stored yet.');
            return;
        end
        [~, roiIdx] = min(abs([roiStruct.frame] - userData.current_frame));
        if ~isempty(userData.roi)
            delete(userData.roi); % Delete the previous ROI if exists
        end

        userData.current_frame = roiStruct(roiIdx).frame;
        set(hSlider, 'Value', userData.current_frame);
        set(hFig, 'UserData', userData);
        update_image();

        userData.roi = drawpolyline('Position', roiStruct(roiIdx).roi); % Restore the closest ROI
        set(hFig, 'UserData', userData);
    end

% Save callback function
    function save_callback(~, ~)
        userData = get(hFig, 'UserData');
        roiStruct = userData.roiStruct;

        outname = strrep(fullFilePath, 'MMStack_Default.ome.tif', 'ROIs.mat');
        save(outname, 'roiStruct');
        disp('ROI structure saved successfully!');

    end


% Update image display
    function update_image()
        userData = get(hFig, 'UserData');
        img = imread(localpath, userData.current_frame);
        img = img(cropPixels:end-cropPixels, cropPixels:end-cropPixels);
        set(hImg, 'CData', img);
        set(hFrameEdit, 'String', num2str(userData.current_frame));
    end
end
