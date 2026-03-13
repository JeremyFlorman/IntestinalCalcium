function [axialSignal,axialBF] = straightenWorm(filepath)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% settings
filepath = 'Z:\OAS\foodEncounter\zfis178_FoodEncounter\tph-1-noFood-0minCycleStart\260305_zfis178_tph-1-noFood-0min_1';
isremote = 1;
startIndex = 1; % which video to start analysis.
startframe = 1; % when to begin analysis

crop = 0; % num pixels to crop from image edge. set to 0 for no cropping.
useautothreshold = 1;% set to 1 to calculate a threshold for each image.
useadaptivethreshold = 0; % if useautothreshold is set to 0 adaptive thresholding can be used
removevignette = 0; % if not zero, size of kernel to use for flatfield correction.

minwormarea = 10000; %lower limit to worm area
maxwormarea = 20000; % upper limit to worm area
numSegments = 400; % number of segments to sample when measuring axial signal
axSigLen = 200; % how many pixels to use for registering axial signal.
axSigHeight = 5; % how many pixels to sample across the width of the worm (i.e. dorsal to ventral)

%%
%%
imgDir = dir([filepath '\**\*behavior\*.h5']);
imgDir = unique({imgDir.folder});


straightenedWormPath = 'C:/users/jeremy/desktop/straightenedWorm.tif';

if exist(straightenedWormPath, 'file')
    delete(straightenedWormPath)
end

for nf =startIndex:length(imgDir)

    path = imgDir{nf}

    [fold, nm, ~] = fileparts(path);
    protopath = regexp(fold,'\', 'split');

    if ~isempty(regexp(protopath{end}, '\w+_\d{1,2}$', 'once')) % does the filename already contain an experiment number?
        expSuffix = protopath{end};
    else
        expSuffix = [protopath{end} '_' num2str(nf)];
    end

    protosavename = [fold '\' expSuffix];


    behIndPath = ['C:\tmp\' nm 'behInd.mat'];
    gcampIndPath = ['C:\tmp\' nm 'gcampInd.mat'];
    logEventPath = ['C:\tmp\' nm 'logEvents.mat'];

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
            load(logEventPath)
        end
    else
        disp('Synching Camera feeds...')
        [behavior_indices, gcamp_indices] = syncDualCameras(path,isremote);
    end


    %% get info from log files
    timestamps = behavior_indices.timestamps;

    if isremote == 1
        if ~isfile(logEventPath)
            disp('Processing Log Files...')
            log_events = processLogFile(path,timestamps, isremote);
            save(logEventPath, 'log_events')
        elseif isfile(logEventPath)
            disp('loading event log')
            load(logEventPath)
        end
    else
        disp('Processing Log Files...')
        log_events = processLogFile(path,timestamps, isremote);
    end

    stimTimes = log_events.stimTimes;
    velocity = log_events.velocity;
    xLoc = log_events.xLoc;
    yLoc = log_events.yLoc;


    %% Image parameters
    imgWidth = behavior_indices.img_size{1}(1);
    imgHeight = behavior_indices.img_size{1}(2);
    nFrames = length(timestamps);

    if imgHeight <500
        axSigHeight = 7;
        SEsize = 5;
        szFilter = 100;
    elseif imgHeight <1000
        axSigHeight = 15;
        SEsize = 8;
        szFilter = 100;
    elseif imgHeight >1000
        axSigHeight = 35;
        SEsize = 20;
        szFilter = 2000;
    end


    SEclose = strel('diamond',SEsize);
    SEopen = strel('diamond', 4);

    axialSignal = NaN(nFrames, round(numSegments/3));
    axialBF = NaN(nFrames, round(numSegments/3));
    bulkSignal = NaN(nFrames,1);
    backgroundSignal = NaN(nFrames,1);
    % background1Pct = NaN(nFrames,1);
    orientation = NaN(nFrames,1);
    area = NaN(nFrames,1);
    wormLength = NaN(nFrames, 1);
    prevTrace = [];
    % axialBrightField = NaN(axSigHeight, axSigLen,nFrames);

    %% Tracking Block
    tic
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

        mCh = behavior_h5(:,:, beh_relative_index);
        GFP = gcamp_h5(:,:,gcamp_relative_index);

        worm = segWorm(mCh,1,0,0);
        %
        % if removevignette ~=0
        %     mCh = imflatfield(mCh,removevignette);
        % end
        %
        % if crop ~= 0
        %     mCh = mCh(crop:end-crop,crop:end-crop);
        %     GFP = GFP(crop:end-crop,crop:end-crop);
        %     imgWidth = size(mCh, 2);
        %     imgHeight = size(mCh,1);
        % end
        % %         imshow(imadjust(mCh))
        % % set up thresholding and binary image processing
        % if useautothreshold ==1
        %     T = graythresh(mCh);
        %     BW = imbinarize(mCh, T); % create mask
        % elseif useadaptivethreshold == 1
        %     BW = imbinarize(mCh,'adaptive','ForegroundPolarity','dark','Sensitivity',0.6);
        % end
        %
        %
        % BW = imcomplement(BW);
        % BW = bwmorph(BW,'clean');
        % % BW = bwmorph(BW,'fill');
        % % tempb = BW;
        % BW = bwareaopen(BW, szFilter);
        %
        %
        % BW = imdilate(BW,SEclose);
        % BW = imerode(BW,SEopen);
        % bw1 = BW;
        % BW = imopen(BW, SEclose);
        % % BW = bwmorph(BW,"spur",inf);
        %
        % %
        % %
        % % tempb2 = BW;
        %
        %
        % % identify connected objects
        % CC = bwconncomp(BW);
        % L = labelmatrix(CC);
        % bwprops = regionprops(L, 'Area', 'Centroid','Orientation');
        %
        % % Filter connected components to get biggest most central object.
        % if ~isempty(bwprops)
        %     % xy = vertcat(bwprops.Centroid);
        %     % x = xy(:,1);
        %     % y = xy(:,2);
        %     % distances = sqrt((imgWidth/2 - x) .^ 2 + (imgHeight/2 - y) .^ 2);
        %     % [~, centralIdx] = min(distances); % most central object
        %
        %     % Precompute once per frame
        %     cx = imgWidth * 0.5;
        %     cy = imgHeight * 0.5;
        %
        %     C = reshape([bwprops.Centroid], 2, []).';
        %     dx = C(:,1) - cx;
        %     dy = C(:,2) - cy;
        %
        %     dist2 = dx.*dx + dy.*dy;   % no sqrt
        %     [~, centralIdx] = min(dist2);% most central object
        %
        %     [~, bigIdx] = max([bwprops.Area]); % largest object
        %
        %     % filtering block: wormIdx is the object that we suspect is the worm.
        %     % if the biggest object is also the most central object than we will
        %     % assume that is the worm. If there is another big object off center,
        %     % as will occur with vignetting, we will check to make sure that the
        %     % most central object is within a size range determined by the values
        %     % of minwormarea and maxwormarea.
        %     if bigIdx == centralIdx && bwprops(bigIdx).Area<= maxwormarea
        %         wormIdx = bigIdx;
        %     elseif bwprops(centralIdx).Area <= maxwormarea && ...
        %             bwprops(centralIdx).Area >= minwormarea
        %         disp(['segmentation error at frame: ' num2str(i)]);
        %         wormIdx = centralIdx;
        %     else
        %         wormIdx = bigIdx;
        %     end
        %
        %     % create a copy of the label matrix Lw that contains only the worm.
        %     Lw = L;
        %
        %     %   imshow(label2rgb(L,'jet','k','shuffle'))
        %     if ~isempty(wormIdx)
        %         Lw(Lw~=wormIdx) = 0;
        %         orientation(i,1) = bwprops(wormIdx).Orientation;
        %
        %         % generate mask, outline and skeleton
        %         mask = logical(Lw);
        %         skel = bwmorph(mask,'thin', inf);
        %         [ep] = bwmorph(skel,'endpoints');


        % Extract Axial signal by sampling perpendicular lines from skeleton %

        % if nnz(ep) >0
        % [ey,ex] = find(ep,1);
        % sortSkel= bwtraceboundary(skel,[ey ex],'E');
        % sortSkel = sortSkel(1:ceil(length(sortSkel)/2),:);

        if ~isempty(worm)
            sortSkel = worm.skeleton.pixels;
            stepSize = 3; % # of points along spine for each spine segment

            % Convert images to double before interpolation
            GFP = double(GFP);
            mCh = double(mCh);

            % Precompute ndgrid
            [Xgrid, Ygrid] = ndgrid(1:size(GFP, 1), 1:size(GFP, 2));

            % Precompute interpolants once
            GFP_interp = griddedInterpolant(Xgrid, Ygrid, GFP, 'linear', 'nearest');
            mCh_interp = griddedInterpolant(Xgrid, Ygrid, mCh, 'linear', 'nearest');


            % Define the desired number of evenly spaced samples
            totalPoints = length(sortSkel);

            % Ensure we don't exceed available points
            segments2Sample = min(numSegments, totalPoints - 1);

            % Evenly select indices along sortSkel
            selectedIndices = round(linspace(1, totalPoints - 1, segments2Sample));

            % Number of points along the perpendicular line
            nPoints = axSigHeight;

            % Initialize matrices for performance
            temptrace = nan(nPoints, segments2Sample);
            tempbf = nan(nPoints, segments2Sample);
            perpX = nan(2, segments2Sample); % Store only two points per segment
            perpY = nan(2, segments2Sample);

            for ii = 1:segments2Sample
                idx = selectedIndices(ii);

                if idx + stepSize < totalPoints
                    seg = sortSkel(idx:idx + stepSize, :);
                else
                    seg = sortSkel(idx:end, :);
                end

                A = [seg(1,2) seg(1,1)]; % [x y] coords of point 1
                B = [seg(end,2) seg(end,1)]; % [x y] coords of point 2

                AB = B - A;
                AB = AB / norm(AB);  % Normalize

                ABperp = AB * [0 -1; 1 0]; % Perpendicular unit vector
                ABmid = (A + B) / 2; % Midpoint

                C = ABmid + axSigHeight * ABperp;
                D = ABmid - axSigHeight * ABperp;

                % Ensure C and D are within bounds
                C = max(min(C, [size(GFP,2), size(GFP,1)]), [1,1]);
                D = max(min(D, [size(GFP,2), size(GFP,1)]), [1,1]);

                perpX(:, ii) = [C(1); D(1)];
                perpY(:, ii) = [C(2); D(2)];

                % Generate linearly spaced points along the perpendicular line
                xq = linspace(C(1), D(1), axSigHeight);
                yq = linspace(C(2), D(2), axSigHeight);


                %Transpose xq and yq for NDGRID format
                temptrace(:, ii) = GFP_interp(yq', xq');
                tempbf(:, ii) = mCh_interp(yq', xq');
            end

            % Convert images back to uint8 (if needed for display later)
            % GFP = uint8(GFP);
            % mCh = uint8(mCh);

            if ~isempty(temptrace)

                tt = resample(max(temptrace), size(axialSignal,2), size(temptrace,2),5,20);  % max?
                abf = resample(mean(tempbf), size(axialSignal,2), size(tempbf,2),5,20);


                % % % % real-time autoFixSignal % % %
                % querryLength = 10; %floor(length(tt)*0.1); % fraction of signal to querry
                % leftMean = mean(abf(1:querryLength),'omitnan');
                % rightMean = mean(abf(length(abf)-querryLength:length(abf)),'omitnan');
                %
                % if leftMean>rightMean
                %     tt = fliplr(tt);
                %     temptrace = fliplr(temptrace);
                %     abf = fliplr(abf);
                %     tempbf = fliplr(tempbf);
                % end



                if i > 1 && ~isempty(prevTrace) % try to use weighted orientation score
                    [totalScore, ~] = scoreOrientation(prevTrace, tt);

                    % threshold / hysteresis
                    if totalScore < -0.5
                        tt = fliplr(tt);
                        temptrace = fliplr(temptrace);
                        abf = fliplr(abf);
                        tempbf = fliplr(tempbf);
                    end

                else % if previous trace unavailable, use regional variance
                    N = length(tt);
                    vf = var(tt(1:floor(0.1*N))); % "Anterior" variance
                    vl = var(tt(floor(0.9*N):N)); % "Posterior" variance

                    varianceScore = vl - vf; % variance score

                    if varianceScore <-50 
                        tt = fliplr(tt);
                        temptrace = fliplr(temptrace);
                        abf = fliplr(abf);
                        tempbf = fliplr(tempbf);
                    end
                end

                    


                % % % % % % % % % % % % % % % % % % % %
                headlength = round(numSegments/3);
                axialBF(i,1:headlength) = mean(tempbf(:,1:headlength),'omitmissing');
                axialSignal(i,1:headlength) = mean(temptrace(:,1:headlength),'omitmissing');



            end

            img = uint8(tempbf(:,1:headlength));
            imwrite(img, 'C:\users\jeremy\desktop\straightedWormHead.tif', 'WriteMode','append')            


            imshow(img, [20 150])
            drawnow

            prevTrace = tt;
        end

        if mod(i,90) == 0
            disp(['Working... ' num2str((i/nFrames)*100) '% complete, just chill...'])
        end
    end

    imwrite(uint8(axialBF)', 'C:/users/jeremy/desktop/testkym.tiff')

    delete(behIndPath)
    delete(gcampIndPath)
    delete(logEventPath)
end