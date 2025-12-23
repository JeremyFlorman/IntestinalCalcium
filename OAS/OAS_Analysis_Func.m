function [] = OAS_Analysis_Func(inputs)
fld = inputs.tiffDir; %'C:\Users\Jeremy\Desktop\220316_zfis178_wildtype_1'; % Folder containing the data you want to analyze

serverfolder = inputs.remoteDir; %'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\itr-1\;'  % upload everything to this location.
if isempty(serverfolder)
    serverfolder = fld;
end

%% settings
startIndex = inputs.startFile; % which video to start analysis.
startframe = inputs.startFrame; % when to begin analysis

uploadresults = inputs.uploadResults; % upload data to remote location (serverfolder)?
isremote = inputs.isRemote;     % is our tiff file on the server? If so, we'll copy to local
% folder to run the analysis then move the results to the
% server.

plotstuff = inputs.plotResults; % display tracking
videostuff =inputs.recordVideo; % record video
framerate = inputs.videoFramerate; % display video/plots every Nth iteration of loop.
fps = inputs.inputFramerate;      % frames per sec of input tiff.
troubleshoot =inputs.troubleshoot; % show binary images instead of regular plots
showNormals = inputs.showNormals;
showAxialSignal = inputs.showAxialSignal;
saveAxialMatrix = 0;

crop = inputs.crop; % num pixels to crop from image edge. set to 0 for no cropping.
useautothreshold = inputs.autoThreshold;% set to 1 to calculate a threshold for each image.
useadaptivethreshold = inputs.adaptiveThreshold; % if useautothreshold is set to 0 adaptive thresholding can be used
removevignette = inputs.flatField; % if not zero, size of kernel to use for flatfield correction.
% loadtiff =inputs.loadTiff; % read entire tiff into memory? faster analysis but requires more ram.


minwormarea = 10000; %lower limit to worm area
maxwormarea = 20000; % upper limit to worm area
numSegments = 100; % number of segments to sample when measuring axial signal
axSigLen = 200; % how many pixels to use for registering axial signal.
axSigHeight = 35; % how many pixels to sample across the width of the worm (i.e. dorsal to ventral)

% SEsize = 20;

%%
%%
imgDir = dir([fld '\**\*behavior\*.h5']);
imgDir = unique({imgDir.folder});


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


    % copy remote files
    % if isremote == 1
    %     logFiles = dir([fld '\*log.txt']);
    % 
    %     for logIdx = 1:length(logFiles)
    %         disp(['Copying log file ' num2str(logIdx) ' of ' num2str(length(logFiles))])
    % 
    %         logPath = fullfile(logFiles(logIdx).folder, logFiles(logIdx).name);
    %         localLogPath = ['C:\tmp\' logFiles(logIdx).name];
    %         copyfile(logPath, localLogPath);
    %     end
    % 
    % 
    %     beh_remote_fp = path;
    %     gcamp_remote_fp = strrep(path, 'behavior', 'gcamp');
    % 
    %     beh_local_fp = strrep(path, fold, 'C:\tmp');
    %     gcamp_local_fp = strrep(beh_local_fp, 'behavior', 'gcamp');
    % 
    %     if ~isfolder(beh_local_fp)
    %         disp('Copying Behavior Files...')
    %         copyfile(beh_remote_fp, beh_local_fp);
    %     end
    % 
    %     if ~isfolder(gcamp_local_fp)
    %         disp('Copying GCaMP Files...')
    %         copyfile(gcamp_remote_fp, gcamp_local_fp);
    %     end
    % 
    %     path = beh_local_fp;
    % end



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

    %% Image parameters
    imgWidth = behavior_indices.img_size{1}(1);
    imgHeight = behavior_indices.img_size{1}(2);
    nFrames = length(timestamps);



    %%
    if plotstuff == 1
        if showAxialSignal == 0
            figure('Position', [978 233 719 653],'Color',[1 1 1]);
            tiledlayout(4,3,'Padding','compact')
            ax1 = nexttile([2 1]);
            ax2 = nexttile([2 1]);
            ax3 = nexttile([2 1]);
            ax4 = nexttile([1 1]);
            ax5 = nexttile([1 1]);
            ax6 = nexttile([1 1]);
            ax7 = nexttile([1 3]);
        elseif showAxialSignal == 1
            figure('Position',[978 233 719 653],'Color',[1 1 1]);
            tiledlayout(9,3,'Padding','compact')
            ax1 = nexttile([3 1]);
            ax2 = nexttile([3 1]);
            ax3 = nexttile([3 1]);
            ax4 = nexttile([2 3]);
            ax7 = nexttile([2 3]);
            velAx = nexttile([2 3]);

        end


        if videostuff == 1
            vidfig = gcf;
            if exist('v','var') == 1
                close(v)
            end
            videopath = [protosavename '_Tracking_Video.mp4'];
            if isremote== 0
                v = VideoWriter(videopath,'MPEG-4');
            else
                localVideoPath = [ 'C:\tmp\' expSuffix '_Tracking_Video.mp4'];
                v = VideoWriter(localVideoPath,'MPEG-4');
            end
            v.FrameRate = 15;
            open(v)
        end
    end

    %     if showNormals ==1
    %     normfig = figure();
    %     normAx = axes(Parent=normfig);
    %     end



    if imgHeight <500
        axSigHeight = 7;
        SEsize = 5;
        szFilter = 100;
    elseif imgHeight <1000
        axSigHeight = 20;
        SEsize = 8;
        szFilter = 100;
    elseif imgHeight >1000
        axSigHeight = 35;
        SEsize = 20;
        szFilter = 2000;
    end


    SEclose = strel('diamond',SEsize);
    SEopen = strel('diamond', 4);

    axialSignal = NaN(nFrames, axSigLen);
    axialBF = NaN(nFrames, axSigLen);
    bulkSignal = NaN(nFrames,1);
    antSignal = NaN(nFrames,1);
    postSignal = NaN(nFrames,1);
    backgroundSignal = NaN(nFrames,1);
    background1Pct = NaN(nFrames,1);
    orientation = NaN(nFrames,1);
    area = NaN(nFrames,1);
    wormLength = NaN(nFrames, 1);

    time = (log_events.time-log_events.time(1))/60; %minutes per frame
    wormIdx = [];

    if saveAxialMatrix == 1
        axialMatrix = NaN(axSigHeight, axSigLen,nFrames);
    end


    hBinary = [];
    hBF = [];
    hGFP = [];
    hKymo = [];
    hBulk = [];
    hBkg = [];
    hAnterior = [];
    hPosterior = [];
    hArea = [];

    hNormals = gobjects(numSegments,1);

    for normIdx = 1:numSegments
        hNormals(normIdx) = line(nan(2,1), nan(2,1),'Color', [0.6 0.6 0.6],'Parent', ax1);
    end

    hAxStim = [];
    hBulkStim = [];
    hAreaStim = [];
    frameError = 0;

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



        if removevignette ~=0
            mCh = imflatfield(mCh,removevignette);
        end

        if crop ~= 0
            mCh = mCh(crop:end-crop,crop:end-crop);
            GFP = GFP(crop:end-crop,crop:end-crop);
            imgWidth = size(mCh, 2);
            imgHeight = size(mCh,1);
        end
        %         imshow(imadjust(mCh))
        % set up thresholding and binary image processing
        if useautothreshold ==1
            T = graythresh(mCh);
            BW = imbinarize(mCh, T); % create mask
        elseif useadaptivethreshold == 1
            BW = imbinarize(mCh,'adaptive','ForegroundPolarity','dark','Sensitivity',0.6);
        end


        BW = imcomplement(BW);
        BW = bwmorph(BW,'clean');
        % BW = bwmorph(BW,'fill');
        tempb = BW;
        BW = bwareaopen(BW, szFilter);


        BW = imdilate(BW,SEclose);
        BW = imerode(BW,SEclose);

        BW = imopen(BW, SEopen);
        % BW = bwmorph(BW,"spur",inf);



        tempb2 = BW;



        if troubleshoot == 1
            imshow(tempb,'Parent', ax2);
            title(ax2,'Initial Threshold');

            imshow(tempb2,'Parent', ax3);
            title(ax3,'Processed Mask');
        end


        % identify connected objects
        CC = bwconncomp(BW);
        L = labelmatrix(CC);
        bwprops = regionprops(L, 'Area', 'Centroid','Orientation');

        % Filter connected components to get biggest most central object.
        if ~isempty(bwprops)
            % xy = vertcat(bwprops.Centroid);
            % x = xy(:,1);
            % y = xy(:,2);
            % distances = sqrt((imgWidth/2 - x) .^ 2 + (imgHeight/2 - y) .^ 2);
            % [~, centralIdx] = min(distances); % most central object

            % Precompute once per frame
            cx = imgWidth * 0.5;
            cy = imgHeight * 0.5;

            C = reshape([bwprops.Centroid], 2, []).';
            dx = C(:,1) - cx;
            dy = C(:,2) - cy;

            dist2 = dx.*dx + dy.*dy;   % no sqrt
            [~, centralIdx] = min(dist2);% most central object

            [~, bigIdx] = max([bwprops.Area]); % largest object

            % filtering block: wormIdx is the object that we suspect is the worm.
            % if the biggest object is also the most central object than we will
            % assume that is the worm. If there is another big object off center,
            % as will occur with vignetting, we will check to make sure that the
            % most central object is within a size range determined by the values
            % of minwormarea and maxwormarea.
            if bigIdx == centralIdx && bwprops(bigIdx).Area<= maxwormarea
                wormIdx = bigIdx;
            elseif bwprops(centralIdx).Area <= maxwormarea && ...
                    bwprops(centralIdx).Area >= minwormarea
                disp(['segmentation error at frame: ' num2str(i)]);
                wormIdx = centralIdx;
            else
                wormIdx = bigIdx;
            end

            % create a copy of the label matrix Lw that contains only the worm.
            Lw = L;


            %   imshow(label2rgb(L,'jet','k','shuffle'))
            if ~isempty(wormIdx)
                Lw(Lw~=wormIdx) = 0;

                orientation(i,1) = bwprops(wormIdx).Orientation;


                % generate mask, outline and skeleton
                mask = logical(Lw);
                outline = bwmorph(mask, 'remove',1);
                %             skel = bwskel(mask,'MinBranchLength', 20);
                skel = bwmorph(mask,'thin', inf);

                outskel = logical(outline+skel);
                [ep] = bwmorph(skel,'endpoints');


                % Extract Axial signal by sampling perpendicular lines from skeleton %

                if nnz(ep) >0
                    [ey,ex] = find(ep,1);
                    sortSkel= bwtraceboundary(skel,[ey ex],'E');
                    sortSkel = sortSkel(1:ceil(length(sortSkel)/2),:);

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

                        try
                            % Fix: Transpose xq and yq for NDGRID format
                            temptrace(:, ii) = GFP_interp(yq', xq');
                            tempbf(:, ii) = mCh_interp(yq', xq');
                        catch
                        end
                    end

                    hold off


                    % Convert images back to uint8 (if needed for display later)
                    GFP = uint8(GFP);
                    mCh = uint8(mCh);
                    hold off






                    if ~isempty(temptrace)
                        rsAxMat = resample(temptrace,axSigLen,size(temptrace,2));

                        tt = resample(max(temptrace), size(axialSignal,2), size(temptrace,2),5,20);  % max?
                        abf = resample(mean(tempbf), size(axialSignal,2), size(tempbf,2),5,20);

                        % % % % real-time autoFixSignal % % %
                        querryLength = length(tt)*0.1; % fraction of signal to querry
                        leftMean = mean(tt(1:querryLength),'omitnan');
                        rightMean = mean(tt(length(tt)-querryLength:length(tt)),'omitnan');

                        if leftMean>rightMean
                            tt = fliplr(tt);
                            temptrace = fliplr(temptrace);
                            abf = fliplr(abf);
                            tempbf = fliplr(tempbf);
                        end

                        % % % % % % % % % % % % % % % % % % % %

                        axialBF(i,1:size(abf,2)) = abf;
                        axialSignal(i,1:size(tt,2)) = tt;
                    end
                end

                %% Bulk signal
                blksig = GFP(mask);
                bulkSignal(i,1) = mean(blksig,"all",'omitnan');

                %% Background 1% - lowest 1% of values outside the ROI
                bkgsig = GFP(~mask);
                thresh = prctile(bkgsig, 1);
                bkgMask = false(size(mask));
                bkgMask(~mask) = GFP(~mask)<=thresh;
                background1Pct(i,1) = mean(GFP(bkgMask),'all','omitnan');

                %% Background Signal - lowest 5% of values outside the ROI
                bkgsig = GFP(~mask);
                backgroundSignal(i,1) = mean(bkgsig,'all','omitnan');

                %% Worm Area
                area(i,1) = bwprops(wormIdx).Area;

                %% Worm Length 
                wormLength(i) = totalPoints;


                % Upsample temptrace and tempbf to match original sortSkel size
                originalIndices = 1:totalPoints - 1;
                upsampledIndices = linspace(1, totalPoints - 1, segments2Sample);

                try
                    temptrace = interp1(upsampledIndices, temptrace', originalIndices, 'linear', 'extrap')';
                    tempbf = interp1(upsampledIndices, tempbf', originalIndices, 'linear', 'extrap')';
                catch
                end

                % Plot Stuff

                if plotstuff == 1
                    if i == startframe || mod(i - startframe, framerate) == 0

                        if troubleshoot == 0
                            try
                                mdiff = size(mCh,2)-size(tempbf,2);

                                if mdiff>0
                                    mpadTrace = padarray(tempbf,[0, ceil(mdiff/2)],0,'both');
                                    mpadTrace = mpadTrace(:,1:size(mCh, 2));
                                elseif mdiff<0
                                    mpaTrace = tempbf(:,abs(mdiff):size(mCh, 2));
                                end

                                bfAdj = [0 1];
                                gfpAdj = [0 0.3];

                                mpad_Outskel = padarray(outskel, [size(mpadTrace,1),0], 'post');
                                mmergedImage = vertcat(mCh, mpadTrace)-50;
                                mmergedOverlay = imoverlay(mmergedImage, mpad_Outskel, [1 0 0]);

                                gdiff = size(GFP,2)-size(temptrace,2);

                                if gdiff>0
                                    gpadTrace = padarray(temptrace,[0, ceil(gdiff/2)],0,'both');
                                    gpadTrace = gpadTrace(:,1:size(GFP, 2));
                                elseif gdiff<0
                                    gpadTrace = temptrace(:,abs(gdiff):size(GFP, 2));
                                end
                                
                                gpad_Outskel = padarray(outskel, [size(gpadTrace,1),0], 'post');
                                gmergedImage = vertcat(GFP, gpadTrace);
                                gmergedOverlay = imoverlay(imadjust(gmergedImage, gfpAdj), gpad_Outskel, [0 1 0]);

                                % plot brightfield and GCaMP images
                                if i == startframe
                                    hBF = imshow(mmergedOverlay,'Parent', ax2);
                                    title(ax2,'Brightfield');

                                    hGFP = imshow(gmergedOverlay,'Parent', ax3);
                                    colormap(ax3, "turbo")
                                    title(ax3,'GCaMP');
                                end
                            catch
                                frameError = 1;
                                backupmCh = imoverlay(imadjust(mCh, bfAdj), outskel, [1 0 0]);
                                backupGFP = imoverlay(imadjust(GFP,gfpAdj), outskel, [0 1 0]);
                            end

                        elseif troubleshoot == 1

                            imshow(tempb,'Parent', ax2);
                            title(ax2,'Initial Threshold');

                            imshow(tempb2,'Parent', ax3);
                            title(ax3,'Processed Mask');
                        end



                        axsig = smoothdata(axialSignal(1:i,:),1,'gaussian',60)'-median(backgroundSignal(1:i),'omitnan');

                        if i == startframe  % plot for the first time
                            %% Binary Mask
                            hold(ax1, 'on')
                            hBinary = imshow(label2rgb(L,'jet','k','shuffle'),'Parent', ax1);
                            hold(ax1, 'off')
                            %% Normal Vectors
                            if showNormals == 1
                                nSegs = size(perpX, 2);

                                set(hNormals(1:nSegs), {'XData'}, squeeze(num2cell(perpX,1))', ...
                                    {'YData'}, squeeze(num2cell(perpY,1))');

                                nanPairs = repmat({[NaN NaN]}, numSegments-nSegs, 1);
                                set(hNormals(nSegs+1:end), {'XData'}, nanPairs, {'YData'}, nanPairs);
                                title(ax1,'Binary Mask + Normal Vectors');
                                uistack(hNormals, 'top')
                            end


                            %% Bright field image
                            hBF = imshow(mmergedOverlay,'Parent', ax2);

                            title(ax2,'Brightfield');

                            %% GCaMP image
                            hGFP = imshow(gmergedOverlay,'Parent', ax3);

                            title(ax3,'GCaMP');

                            %% Axial Signal
                            hKymo = imagesc(axsig,'Parent',ax4);
                            ax4.CLim = [0 60];
                            ax4.XLim = [1, length(axialSignal)];
                            ax4.XAxis.Visible = 0;
                            ax4.YTickLabel = [];
                            ax4.YTick = [];
                            % ylabel(ax4, 'Longitudinal Kymograph')
                            box(ax4, 'off')
                            colormap turbo
                            cb = colorbar(ax4);
                            cb.Location = 'manual';
                            cb.Position = [0.0567 0.4828 0.0064 0.1594];
                            cb.Label.String = 'Fluorescent Intensity (a.u.)';


                            %% Bulk signal
                            hBulk = plot(time,bulkSignal, 'Parent', ax7);
                            hold(ax7, 'on')
                            hBkg = plot(time',backgroundSignal, 'Parent', ax7);
                            hAnterior = plot(time, antSignal,'g', 'Parent', ax7);
                            hPosterior = plot(time, postSignal, 'r', 'Parent', ax7);

                            hold(ax7, 'off')
                            ax7.XLim = [0 time(end)];
                            ylabel(ax7, 'Bulk Ca^2^+ Signal');
                            box(ax7, 'off')
                            ax7.TickLength = [0.005 0.005];
                            
                            %% Area
                            hArea = plot(time(1:i),smoothdata(area(1:i),'gaussian', 30), 'Parent', velAx);
                            xlim([0 time(end)]);
                            ylabel(velAx, 'Worm Area (Pixels)');
                            xlabel(velAx,'Time (min)');
                            velAx.TickLength = [0.005 0.005];
                            box off
                        end

                        %% plot stimuli
                        if ~isempty(stimTimes)
                            stimX = stimTimes(stimTimes<=i);
                            axStimY = repmat(1, size(stimX));
                            areaStimY = repmat(velAx.YLim(2)*0.99, size(stimX));
                            bulkStimY = repmat(ax7.YLim(2)*0.99, size(stimX));

                            if i >= stimTimes(1) && i <= stimTimes(1)+framerate
                                hAxStim = line(stimX,axStimY,'Marker', 'diamond', 'Marker', 'v', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0],'LineStyle', 'none','Parent', ax4);
                                hBulkStim = line(time(stimX),bulkStimY,'Marker', 'v', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0], 'LineStyle', 'none', 'Parent', ax7);
                                hAreaStim = line(time(stimX),areaStimY,'Marker', 'v', 'MarkerSize', 8, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0],'LineStyle', 'none','Parent', velAx);
                            end
                        end

                        %% update plots
                        if i ~= startframe
                            % Binary Mask
                            hBinary.CData = label2rgb(L,'jet','k','shuffle');
                            numSegs = size(perpX,2);
                            set(hNormals(1:numSegs), {'XData'}, squeeze(num2cell(perpX,1))', ...
                                {'YData'}, squeeze(num2cell(perpY,1))');

                            nanPairs = repmat({[NaN NaN]}, numel(hNormals)-numSegs, 1);
                            set(hNormals(numSegs+1:end), {'XData'}, nanPairs, {'YData'}, nanPairs);

                            if troubleshoot == 0
                                if frameError == 0
                                % Brightfield
                                hBF.CData = mmergedOverlay;

                                % GCaMP
                                hGFP.CData = gmergedOverlay;
                                elseif frameError == 1
                                    hBF.CData = backupmCh;
                                    hGFP.CData = backupGFP;
                                    frameError = 0;
                                end
                            end

                            % Axial Signal
                            hKymo.CData = axsig;

                            % Bulk & Background Signal
                            hBulk.XData = time;
                            hBulk.YData = bulkSignal;

                            hAnterior.XData = time;
                            hAnterior.YData = antSignal;

                            hPosterior.XData = time;
                            hPosterior.YData = postSignal;

                            hBkg.XData = time';
                            hBkg.YData = backgroundSignal;

                            % Area
                            hArea.XData = time(1:i);
                            hArea.YData = smoothdata(area(1:i),'gaussian', 30);

                            if ~isempty(stimTimes) && i>stimTimes(1)
                                hAxStim.XData = stimX;
                                hAxStim.YData = axStimY;

                                hBulkStim.XData = time(stimX);
                                hBulkStim.YData = bulkStimY;

                                hAreaStim.XData = time(stimX);
                                hAreaStim.YData = areaStimY;
                            end


                            drawnow limitrate
                        end

                        if videostuff == 1
                            frame = getframe(vidfig);
                            writeVideo(v,frame);
                        end
                    end
                end
            end
        end
        if mod(i,90) == 0
            disp(['Working... ' num2str((i/nFrames)*100) '% complete, just chill...'])
            % m = memory;
            % disp(['Memory Usage: ' num2str(m.MemUsedMATLAB/1073741824) ' Gb'])
        end
    end

    disp('file processed in:')
    toc

    if exist('v','var') == 1
        close(v)
    end
    %% %%%%%%%%%%%%%%%%%%%%% Auto Fix Axial Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fractionToQuerry = 0.1;
    % autoAxialSignal = autoFixSignal(axialSignal,fractionToQuerry);
    autoAxialSignal = axialSignal;

    %     for ii = 1:length(axialSignal)
    %         left = mean(axialSignal(ii,1:10),'omitnan');
    %         right = mean(axialSignal(ii,end-10:end),'omitnan');
    %         if left > right
    %             autoAxialSignal(ii,:) = fliplr(axialSignal(ii,:));
    %             autoAxialBF(ii,:) = fliplr(axialBF(ii,:));
    %             axmat(ii,1) = {fliplr(axmat{ii,1})};
    %             axmat(ii,2) = {fliplr(axmat{ii,2})};
    %         elseif left <= right
    %             autoAxialSignal(ii,:) = axialSignal(ii,:);
    %             autoAxialBF(ii,:) = axialBF(ii,:);
    %         end
    %     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    % peak analysis
    [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',3,'MinPeakDistance',150);
    peakpad = fps*15; % framerate*time in seconds;
    pktime = linspace(-15,15, peakpad*2)';
    if isempty(loc)
        loc = NaN;
        locinmin = NaN;
    end

    pktraces = NaN(peakpad*2,length(loc));
    for k = 1:length(loc)
        prepeak = loc(k)-peakpad;
        postpeak = loc(k)+peakpad-1;
        temppeak = [];

        if prepeak>1 && postpeak<length(bulkSignal)
            pktraces(1:peakpad*2,k) = bulkSignal(prepeak:postpeak,1);
        elseif prepeak>1 && postpeak>length(bulkSignal)
            temptrace = bulkSignal(prepeak:end,1);
            pktraces(1:length(temptrace),k) = temptrace;
        end

    end
    pkmean = mean(pktraces,2,'omitnan');

    %% Save Stuff

    datasavename = [protosavename '_wormdata.mat']


    wormdata = struct();
    wormdata.autoAxialSignal = autoAxialSignal;
    % wormdata.sumSignal = sumSignal;
    wormdata.bulkSignal = bulkSignal;
    % wormdata.bulkAboveBkg = bulkAboveBkg;
    wormdata.backgroundSignal = backgroundSignal;
    wormdata.background1Pct = background1Pct;
    wormdata.orientation = orientation;
    wormdata.area = area;
    wormdata.wormLength = wormLength;
    wormdata.peakTraces = pktraces;
    wormdata.peakLoc = loc;
    wormdata.include = 1;
    wormdata.stimTimes = stimTimes;
    wormdata.velocity = log_events.velocity;

    save(datasavename, 'wormdata')

    %% load stuff
    % tic
    %     load(datasavename)
    %     autoAxialSignal = wormdata.autoAxialSignal;
    %     sumSignal = wormdata.sumSignal;
    %     bulkSignal = wormdata.bulkSignal;
    %     bulkAboveBkg = wormdata.bulkAboveBkg;
    %     backgroundSignal =wormdata.backgroundSignal;
    %     wormdata.orientation = orientation;
    %     area = wormdata.area;
    %     pktraces = wormdata.peakTraces;
    %     loc = wormdata.peakLoc;
    %     stimTimes = wormdata.stimTimes;
    %     velocity = wormdata.velocity;

    %% Plot traces
    if ~exist('time','var')
        time = linspace(0,round((nFrames)/fps/60,1),nFrames); %minutes per frame
    end
    if ~exist('pk','var')
        [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',3, 'MinPeakDistance',150);
        peakpad = fps*15;
        pktime = linspace(-15,15, peakpad*2)';
        pkmean = mean(pktraces,2,'omitnan');
    end



    figure('Position', [936 72 903 586],Color=[1 1 1])
    t = tiledlayout(4,4,'TileSpacing','compact','Padding','tight');

    % % % Bulk Signal % % %
    nexttile([1 3])
    if ~isnan(loc)
        plot(time,bulkSignal-backgroundSignal,time(loc),pk*1.01, 'rv')
    else
        plot(time,bulkSignal-backgroundSignal)
    end
    hold on
    if ~isempty(stimTimes)
        ax = gca;
        plot(time(stimTimes),ax.YLim(2)*.98,'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
    end

    % plot(time, backgroundSignal)
    hold off

    xlim([0 time(end)])
    %     xlabel(gca, 'Time (min)')
    ylabel(gca,'Fluorescence (a.u.)');
    title(gca, 'Whole Animal Calcium Trace')
    ax = gca;
    xt = ax.XTick;
    xtl = ax.XTickLabels;
    ax.TickLength =[0.005 0.005];
    box off

    % % % Peak Profile % % %
    nexttile([1 1])
    plot(pktime, pktraces, 'Color', [0.7 0.7 0.7])
    hold on
    plot(pktime, pkmean, 'Color', [1 0 0], 'LineWidth', 2);
    hold off
    title(gca, 'Spike Profile');
    ylabel(gca,'Fluorescence (a.u.)');
    colormap bone
    box off

    % % % Axial Signal % % %
    ax = nexttile([1 3]);
    imagesc(smoothdata(autoAxialSignal,1,'gaussian',60)'-median(backgroundSignal,'omitnan'))
    title(gca, 'Axial Calcium Trace')
    hold on
    plot(loc,1, 'vw', 'MarkerFaceColor' ,[.4 .5 .6]);
    if ~isempty(stimTimes)
        plot(stimTimes,1,'Marker', 'diamond', 'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
    end
    hold off
    box off

    %     xlabel('Time (min)')
    ax.XTick = xt*60*fps; %linspace(0,length(autoAxialSignal),length(xtl));
    ax.XTickLabels = xtl;
    ax.YTick = [20 size(autoAxialSignal,2)-20];
    ax.YTickLabel = {'Head', 'Tail'};
    ax.CLim =[0 50];
    colormap turbo
    ax.TickLength = [0.001 0.001];


    % % % Interval Histogram % % %
    nexttile([1 1])
    edges = 0:2:120;
    histogram(diff(loc)./fps,'BinEdges',edges);
    title(gca,'Inter-Peak Interval');
    % ylim([0 10])
    xlim([0 120])
    xlabel(gca,'Time (s)');
    ylabel(gca,'Count');
    box off

    % % % Velocity  % % %
    nexttile([1 3]);
    plot(time,smoothdata(velocity,'gaussian',30))
    ax = gca;
    hold on
    if ~isempty(stimTimes)
        plot(time(stimTimes),ax.YLim(2)*.98,'Marker', 'v', 'MarkerSize', 9, 'MarkerFaceColor', [0.8 .2 .5], 'MarkerEdgeColor', [0 0 0])
    end
    hold off
    xlim([0 time(end)])
    title(gca, 'Velocity')
    ylabel(gca,'Steps/sec')
    %     xlabel(gca,'Time (min)')
    ax.TickLength = [0.005 0.005];
    box off

    % % % Peak Widths % % %
    nexttile([1 1])
    histogram(w./fps,'BinEdges', 1:15);
    ylim([0 10])
    xlim([0 15])
    title(gca,'Peak Widths');
    ylabel(gca,'Count');
    xlabel(gca,'Time (s)');
    box off
    % % % Worm Area % % %
    nexttile([1 3]);
    plot(time,smoothdata(area,'gaussian', 30))
    xlim([0 time(end)])
    title(gca, 'Worm Area')
    ylabel(gca,'Pixels')
    xlabel(gca,'Time (min)')
    ax =  gca;
    ax.TickLength = [0.005 0.005];
    box off


    %     nexttile([1 1]);
    %     plot(orientation,time);
    %     ylim([0 time(end)])
    %     title(gca, 'Worm Orientation')
    %     set(gca,'YDir', 'reverse')
    %     xlabel(gca,'Orientation (degrees)')
    %     ylabel(gca, 'Time (min)')


    title(t, strrep(expSuffix, '_', ' '));

    summaryPlotName = [protosavename '_Summary_Plots.png'];

    saveas(gcf, summaryPlotName)



    %% Copy to server
    if uploadresults == 1
        if isremote == 0  % if working with local files, upload to serverfolder (specified in settings)
            [folder2Copy, ~]= fileparts(path);
            lastfolder = regexp(folder2Copy, '\', 'split');
            serverLocation = [serverfolder '\' expSuffix];

            if ~isfolder(serverLocation)
                mkdir(serverLocation);
            end


            % copy behavior h5 files
            behaviorH5Dir = imgDir{nf};
            [parentfolder, h5folder] = fileparts(behaviorH5Dir);
            [statusbeh,~,~]=copyfile(behaviorH5Dir, [serverLocation '\' h5folder '\']);

            % copy GCaMP h5 files
            gcampH5Dir = strrep(behaviorH5Dir, 'behavior', 'gcamp');
            [statusgc,~,~]=copyfile(gcampH5Dir, [serverLocation '\' strrep(h5folder, 'behavior', 'gcamp') '\']);

            % copy log files

            logDir = dir([parentfolder '\*log.txt']);
            for li = 1:length(logDir)
                [statuslog,~,~]=copyfile(fullfile(logDir(li).folder,logDir(li).name), serverLocation);
            end

            % copy wormdata
            [statuswormdata,~,~]=copyfile(datasavename, serverLocation);

            % copy summary plots
            [statussummaryplot,~,~]=copyfile(summaryPlotName, serverLocation);

            % copy summary plots
            [statusvideoplot,~,~]=copyfile(videopath, serverLocation);


        elseif isremote == 1  % if working with remote files, moved analyzed results back to where we found them.
            clear('img')
            [statusvideoplot,~,~]=copyfile(localVideoPath, videopath);
            delete(localVideoPath);
        end
    end

    if exist('wormdata', 'var')
        clear('wormdata');
    end

    clear('h5Data')
    clear('gfp')
    clear('bf')


    if exist('img', 'var')
        clear('img')
    end

    if isremote == 1
        delete(behIndPath)
        delete(gcampIndPath)
        delete(logEventPath)
    end

end
