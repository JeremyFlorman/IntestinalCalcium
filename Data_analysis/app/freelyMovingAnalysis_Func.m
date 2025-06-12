function [] = freelyMovingAnalysis_Func(inputs)
fld = inputs.tiffDir; %'C:\Users\Jeremy\Desktop\220316_zfis178_wildtype_1'; % Folder containing the data you want to analyze
serverfolder = inputs.remoteDir; %'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\itr-1\;'  % upload everything to this location.

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
loadtiff =inputs.loadTiff; % read entire tiff into memory? faster analysis but requires more ram.


minwormarea = 1000; %lower limit to worm area
maxwormarea = 4000; % upper limit to worm area
numSegments = 100; % number of segments to sample when measuring axial signal
axSigLen = 200; % how many pixels to use for registering axial signal.(i.e. pixels from head to tail)
axSigHeight = 10; % how many pixels to sample across the width of the worm (i.e. dorsal to ventral)
SEsize = 20;


%%
tdir = dir([fld '\**\*.tif']);

for nf =startIndex:length(tdir)
    path = fullfile(tdir(nf).folder, tdir(nf).name)

    if isremote == 1

        uploadresults = 1;
        tempfolder = 'C:\tmp'; %set temp directory for copying tiff files.


        if ~exist(fullfile(tempfolder, tdir(nf).name), "file") % copy tiff if it isn't already in tempfolder;
            tic
            disp(['Copying file to: ' tempfolder ' for local processing'])
            copyfile(path, tempfolder)
            disp(['File copied in ' num2str(toc)])
        end

        remotepath = path;                         % save the original path for later uploading.
        path = fullfile(tempfolder, tdir(nf).name); % rename path to the local path.

    end

    %%
    if plotstuff == 1
        if showAxialSignal == 0
            figure('Position', [668 128 1064 653],'Color',[1 1 1]);
            tiledlayout(4,3,'Padding','compact')
            ax1 = nexttile([2 1]);
            ax2 = nexttile([2 1]);
            ax3 = nexttile([2 1]);
            ax4 = nexttile([1 1]);
            ax5 = nexttile([1 1]);
            ax6 = nexttile([1 1]);
            ax7 = nexttile([1 3]);
        elseif showAxialSignal == 1
            figure('Position',[668 128 1064 653],'Color',[1 1 1]);
            tiledlayout(4,3,'Padding','compact')
            ax1 = nexttile([2 1]);
            ax2 = nexttile([2 1]);
            ax3 = nexttile([2 1]);
            ax4 = nexttile([1 3]);
            ax7 = nexttile([1 3]);
        end


        if videostuff == 1
            vidfig = gcf;
            if exist('v','var') == 1
                close(v)
            end
            videopath = strrep(path, '_MMStack_Default.ome.tif', '_Tracking_Video.mp4');
            v = VideoWriter(videopath,'MPEG-4');
            v.FrameRate = 15;
            open(v)
        end
    end

    tic
    info = imfinfo(path);
    toc
    imgWidth = info(1).Width;
    imgHeight = info(1).Height;



    axialSignal = NaN(length(info)/2, axSigLen);

    axialBF = NaN(length(info)/2, axSigLen);
    % sumSignal = NaN(length(info)/2,1);
    bulkSignal = NaN(length(info)/2,1);
    % bulkAboveBkg = NaN(length(info)/2,1);
    backgroundSignal = NaN(length(info)/2,1);
    orientation = NaN(length(info)/2,1);
    area = NaN(length(info)/2,1);
    wormLength = NaN(length(info)/2,1);


    time = linspace(0,round((length(info)/2)/fps/60,1),ceil(length(info)/2)); %minutes per frame
    wormIdx = [];

    if saveAxialMatrix == 1
        axialMatrix = NaN(axSigHeight+1, axSigLen,length(info)/2);
    end


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    roiPath = dir([tdir(nf).folder '\*ROIs.mat']);
    if ~isempty(roiPath)
        useROI = 1;
        load(fullfile(roiPath.folder,roiPath.name));
        % roiFrame = regexpi(roiPath.name, 'roi_(\d+).mat','tokens');
        % roiFrame = str2double(roiFrame{1});
        roiFrames = ceil([roiStruct.frame]/2);
        imageSize = [251, 251];  % Replace with the size of your image
        binaryMask = false(imageSize(1), imageSize(2), length(roiFrames));

        for i = 1:length(roiFrames)
            if ~isempty(roiStruct(i).roi)
                x = roiStruct(i).roi(:,1);
                y = roiStruct(i).roi(:,2);
                numPoints = 200;
                % Interpolate points along the polyline
                xInterp = interp1(1:numel(x), x, linspace(1, numel(x), numPoints), 'linear');
                yInterp = interp1(1:numel(y), y, linspace(1, numel(y), numPoints), 'linear');

                % Round to the nearest pixel coordinates
                xInterp = round(xInterp);
                yInterp = round(yInterp);


                bm = false(imageSize);

                % Set the interpolated pixel coordinates in the binary mask to 1
                for k = 1:length(xInterp)
                    if xInterp(k) >= 1 && xInterp(k) <= imageSize(2) && yInterp(k) >= 1 && yInterp(k) <= imageSize(1)
                        bm(yInterp(k), xInterp(k)) = true;
                    end
                end

                bm = imdilate(bm,strel('disk',2));
                binaryMask(1:imageSize(1), 1:imageSize(2),i) = bwmorph(bm, 'thin', inf);
            end
        end

    else
        useROI = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    if loadtiff == 1
        if exist('TIFFStack\','dir')
            img = TIFFStack(path);
            war = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
            warning('off', 'MATLAB:imagesci:tifftagsread:expectedTagDataFormat');
            warning('off','imageio:tiffmexutils:libtiffWarning')
            warning('off','imageio:tiffutils:libtiffWarning')
        else
            disp('This Code works way faster with the TiffStack function: https://github.com/DylanMuir/TIFFStack')
            img = tiffreadVolume(path);
        end
        disp('loading tif took:')
        toc
    end

    %% Tracking Block
    for i = startframe:length(info)/2
        if loadtiff == 1
            mCh = img(:,:,(i*2)-1);      % mCherry signal is odd # frames?
            GFP = img(:,:, i*2);    % GCaMP signal is even # frames?
        elseif loadtiff == 0
            mCh = imread(path, (i*2)-1);      % mCherry signal is odd # frames?
            GFP = imread(path, i*2);    % GCaMP signal is even # frames?
        end


        if removevignette ~=0
            mCh = imflatfield(mCh,removevignette);
        end

        if crop ~= 0
            mCh = mCh(crop:end-crop,crop:end-crop);
            GFP = GFP(crop:end-crop,crop:end-crop);
            imgWidth = size(mCh, 2);
            imgHeight = size(mCh,1);
        end

        % set up thresholding and binary image processing
        if useautothreshold ==1
            T = graythresh(mCh);
            BW = imbinarize(mCh, T); % create mask
        elseif useadaptivethreshold == 1
            BW = imbinarize(mCh,'adaptive','ForegroundPolarity','dark','Sensitivity',0.6);
        end


        BW = imcomplement(BW);
        % BW = bwmorph(BW,'clean');
        % BW = bwmorph(BW,'fill');
        tempb = BW;

        
        BW = bwareaopen(BW, 100);
        BW = imfill(BW,'holes');
        
        BW = imdilate(BW,strel('disk',6));
        BW = imerode(BW,strel('disk',6));
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
            xy = vertcat(bwprops.Centroid);
            x = xy(:,1);
            y = xy(:,2);
            distances = sqrt((imgWidth/2 - x) .^ 2 + (imgHeight/2 - y) .^ 2);
            [~, centralIdx] = min(distances); % most central object
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

                if useROI == 1
                    roiIndex = find(roiFrames<= i,1,'last');
                end
                
                if useROI == 1 && i>= roiFrames(1) && nnz(binaryMask(:,:,roiIndex))>0
                    skel = binaryMask(:,:,roiIndex);
                    roiActive = 1;
                else
                    skel = bwmorph(mask,'thin', inf);
                    roiActive = 0;
                end

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
                    wormLength(i) = totalPoints;

                    % Ensure we don't exceed available points
                    segments2Sample = min(numSegments, totalPoints - 1);

                    % Evenly select indices along sortSkel
                    selectedIndices = round(linspace(1, totalPoints - 1, segments2Sample));

                    % Number of points along the perpendicular line
                    if roiActive == 0
                        nPoints = axSigHeight;
                    elseif roiActive == 1
                        nPoints = 4;
                    end

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
                    GFP = uint16(GFP);
                    mCh = uint16(mCh);
                    hold off




                    % temptrace = cell2mat(temptrace);
                    % tempbf = cell2mat(tempbf);
                    % perpX = cell2mat(perpX);
                    % perpY = cell2mat(perpY);


                    if ~isempty(temptrace)
                        try
                            % resample axial images
                            x1 = linspace(1,size(temptrace,2), size(temptrace,2));
                            x2 = linspace(1,size(temptrace,2),axSigLen);
                            temptrace = interp1(x1, temptrace',x2)';
                            tempbf = interp1(x1, tempbf',x2)';

                            tt = max(temptrace);
                            abf = mean(tempbf);


                            % % % % real-time autoFixSignal % % %
                            % if roiActive == 0
                                querryLength = length(tt)*0.1; % fraction of signal to querry
                                leftMean = mean(tt(1:querryLength),'omitnan');
                                rightMean = mean(tt(length(tt)-querryLength:length(tt)),'omitnan');

                                if leftMean>rightMean
                                    tt = fliplr(tt);
                                    temptrace = fliplr(temptrace);
                                    abf = fliplr(abf);
                                    tempbf = fliplr(tempbf);
                                end
                            % end
                            % % % % % % % % % % % % % % % % % % % %

                            axialBF(i,1:size(abf,2)) = abf;
                            axialSignal(i,1:size(tt,2)) = tt;

                            if saveAxialMatrix == 1
                                axialMatrix(:,:,i) = temptrace;
                            end
                        catch
                        end
                    end
                end

                % Bulk signal and background signal
                blksig = GFP(mask);
                % sumSignal(i,1) = sum(blksig,"all",'omitnan');
                bulkSignal(i,1) = mean(blksig,"all",'omitnan');
                backgroundSignal(i,1) = mean(GFP(~mask),'all','omitnan');

                % abovebkg = blksig>mean(GFP(~mask));
                % bulkAboveBkg(i,1) = mean(blksig(abovebkg)-mean(GFP(~mask)),"all",'omitnan');

                area(i,1) = bwprops(wormIdx).Area;



                if plotstuff == 1
                    if mod(i,framerate) == 0

                        imshow(label2rgb(L,'jet','k','shuffle'),'Parent', ax1)
                        
                        if showNormals == 1
                            line(perpX,perpY,'Color', [0.9 0.9 0.9],'Parent', ax1)
                            title(ax1,'Binary Mask + Normal Vectors');
                        else 
                            title(ax1,'Binary Mask');
                        end

                        if troubleshoot == 1
                            imshow(tempb,'Parent', ax2);
                            title(ax2,'Initial Threshold');

                            imshow(tempb2,'Parent', ax3);
                            title(ax3,'Processed Mask');
                        else
                            try

                                mdiff = size(mCh,2)-size(tempbf,2);

                                if mdiff>0
                                    mpadTrace = padarray(tempbf,[0, ceil(mdiff/2)],0,'both');
                                    mpadTrace = mpadTrace(:,1:size(mCh, 2));
                                elseif mdiff<0
                                    mpaTrace = tempbf(:,abs(mdiff):size(mCh, 2));
                                end

                                mpad_Outskel = padarray(outskel, [size(mpadTrace,1),0], 'post');
                                mmergedImage = vertcat(mCh, mpadTrace);
                                mmergedOverlay = imoverlay(imadjust(mmergedImage, [0.05 0.21]), mpad_Outskel, [1 0 0]);
                                imshow(mmergedOverlay,'Parent', ax2)
                                title(ax2,'Brightfield');

                                gdiff = size(GFP,2)-size(temptrace,2);

                                if gdiff>0
                                    gpadTrace = padarray(temptrace,[0, ceil(gdiff/2)],0,'both');
                                    gpadTrace = gpadTrace(:,1:size(GFP, 2));
                                elseif gdiff<0
                                    gpadTrace = temptrace(:,abs(gdiff):size(GFP, 2));
                                end

                                gpad_Outskel = padarray(outskel, [size(gpadTrace,1),0], 'post');
                                gmergedImage = vertcat(GFP, gpadTrace);
                                gmergedOverlay = imoverlay(imadjust(gmergedImage, [0.06 0.2]), gpad_Outskel, [0 1 0]);
                                imshow(gmergedOverlay,'Parent', ax3)
                                title(ax3,'GCaMP');
                            catch
                                imshow(imoverlay(imadjust(mCh, [0.05 0.21]), outskel, [1 0 0]), 'Parent', ax2)
                                title(ax2,'Brightfield');

                                imshow(imoverlay(imadjust(GFP,[0.06 0.1]), outskel, [0 1 0]), 'Parent', ax3)
                                title(ax3,'GCaMP');
                            end
                        end

                        if showAxialSignal == 0
                            plot(time,area(:), 'Parent', ax4);
                            title(ax4,'Worm Area');
                            ylabel(ax4,'Pixels');
                            xlabel(ax4,'Time (min)');

                            plot(time,orientation(:), 'Parent', ax5);
                            title(ax5,'Worm Orientation');
                            ylabel(ax5,'Degrees');
                            xlabel(ax5,'Time (min)');

                            plot(1:size(axialSignal,2),axialSignal(i,:), 'g',...
                                1:size(axialBF,2),axialBF(i,:),'r', 'Parent', ax6)
                            ax6.XLim = [0 size(axialSignal,2)];
                            ax6.YLim = [0 50000];
                            title(ax6,'Signal Along Midline');
                            ylabel(ax6,'Mean Fluorescent Intensity (a.u.)');
                            xlabel(ax6, 'head <--- Distance (pixels) ---> tail');
                            legend(ax6,{'GCaMP6', 'Brightfield'}, 'Location', 'northwest', ...
                                'Box', 'off');

                        elseif showAxialSignal == 1
                            axsig = smoothdata(axialSignal,1,'gaussian',60)'-median(backgroundSignal(1:i),'omitnan');
                            imagesc(axsig,'Parent',ax4)
                            ax4.CLim = [-500 30000];
                            ax4.XAxis.Visible = 0;
                            ax4.YAxis.Visible = 0;
                            title(ax4,'Longitudinal Kymograph')
                            colormap turbo
                            cb = colorbar(ax4);
                            cb.Location = 'manual';
                            cb.Position = [0.0507 0.3054 0.0095 0.1676];
                            cb.Label.String = 'Fluorescent Intensity (a.u.)';

                        end

                        plot(time,bulkSignal(:),time,backgroundSignal(:), 'Parent', ax7)
                        if i>1
                            xlim(ax7,[0 time(end)]);
                        end
                        title(ax7, 'Mean Bulk Signal')
                        ylabel(ax7, 'Mean Fluorescence (a.u.)');
                        xlabel(ax7,'Time (min)');
                        box off
                        drawnow

                        if videostuff == 1
                            frame = getframe(vidfig);
                            writeVideo(v,frame);
                        end
                    end
                end
            end
            if mod(i,90) == 0
                disp(['Working... ' num2str((i/length(info))*200) '% complete, just chill...'])
            end
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
    [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',750,'MinPeakDistance',150);
    peakpad = fps*15; % framerate*time in seconds;
    pktime = linspace(-15,15, peakpad*2)';
    if isempty(loc)
        loc = NaN;
        locinmin = NaN;
    else
        locinmin = loc/(fps*60);
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
    [folder, name, ext] = fileparts(path);
    outname = fullfile(folder, name);
    datasavename = strrep(outname, 'MMStack_Default.ome', 'wormdata.mat');


    wormdata = struct();
    wormdata.autoAxialSignal = autoAxialSignal;
    if saveAxialMatrix == 1
        wormdata.axialMatrix = axialMatrix;
    end

    % if useROI
    %     wormdata.noAutoFix = 1;
    % end
    % wormdata.sumSignal = sumSignal;
    wormdata.bulkSignal = bulkSignal;
    % wormdata.bulkAboveBkg = bulkAboveBkg;
    wormdata.backgroundSignal = backgroundSignal;
    wormdata.orientation = orientation;
    wormdata.area = area;
    wormdata.wormLength = wormLength;
    wormdata.peakTraces = pktraces;
    wormdata.peakLoc = loc;
    wormdata.include = 1;

    save(datasavename, 'wormdata')

    %% load stuff
    % tic
    % [folder, name, ext] = fileparts(path);
    % outname = fullfile(folder, name);
    % axialSignal= load([outname '_axialSignal.txt']);
    % bulkSignal = load([outname '_bulkSignal.txt']);
    % orientation = load([outname '_orientation.txt']);
    % area=load([outname '_worm_area.txt']);
    % pktraces = load([outname '_peakTraces.txt']);
    % loc = load([outname '_peakLocs.txt']);
    % toc
    %% Plot traces
    if ~exist('time','var')
        time = linspace(0,round((length(info)/2)/fps/60,1),ceil(length(info)/2)); %minutes per frame
    end
    if ~exist('pk','var')
        [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',1000, 'MinPeakDistance',150);
        peakpad = fps*15;
        pktime = linspace(-15,15, peakpad*2)';
        pkmean = mean(pktraces,2,'omitnan');
    end



    figure('Position', [135.4000 142.6000 1232 586.4000],Color=[1 1 1])
    t = tiledlayout(3,4,'TileSpacing','compact','Padding','tight');

    nexttile([1 3])
    if ~isnan(locinmin)
        plot(time,bulkSignal,locinmin,pk*1.01, 'rv')
    else
        plot(time,bulkSignal)
    end

    hold on
    plot(time, backgroundSignal)
    hold off

    xlim([0 time(end)])
    xlabel(gca, 'Time (min)')
    ylabel(gca,'Fluorescence (a.u.)');
    title(gca, 'Whole Animal Calcium Trace')
    ax = gca;
    xtl = ax.XTickLabels;

    nexttile([1 1])
    plot(pktime, pktraces, 'Color', [0.7 0.7 0.7])
    hold on
    plot(pktime, pkmean, 'Color', [1 0 0], 'LineWidth', 2);
    hold off
    title(gca, 'Spike Profile');
    ylabel(gca,'Fluorescence (a.u.)');
    colormap bone


    ax = nexttile([1 3]);
    imagesc(smoothdata(autoAxialSignal,1,'gaussian',60)'-median(backgroundSignal(1:i),'omitnan'))
    title(gca, 'Axial Calcium Trace')
    hold on
    plot(loc,1, 'vw', 'MarkerFaceColor' ,[.4 .5 .6]);
    hold off
    box off

    xlabel('Time (min)')
    ax.XTick = linspace(1,length(autoAxialSignal),length(xtl));
    ax.XTickLabels = xtl;
    ax.YTick = [20 size(autoAxialSignal,2)-20];
    ax.YTickLabel = {'Head', 'Tail'};
    ax.CLim =[-500 30000];
    colormap turbo



    nexttile([1 1])
    edges = 0:2:120;
    histogram(diff(loc)./fps,'BinEdges',edges);
    title(gca,'Inter-Peak Interval');
    ylim([0 10])
    xlim([0 120])
    xlabel(gca,'Time (s)');
    ylabel(gca,'Count');


    nexttile([1 3]);
    plot(time,smoothdata(area,1,'sgolay',fps))
    xlim([0 time(end)])
    title(gca, 'Worm Area')
    ylabel(gca,'Pixels')
    xlabel(gca,'Time (min)')

    nexttile([1 1])

    histogram(w./fps,'BinEdges', 1:15);
    ylim([0 10])
    xlim([0 15])
    title(gca,'Peak Widths');
    ylabel(gca,'Count');
    xlabel(gca,'Time (s)');



    %     nexttile([1 1]);
    %     plot(orientation,time);
    %     ylim([0 time(end)])
    %     title(gca, 'Worm Orientation')
    %     set(gca,'YDir', 'reverse')
    %     xlabel(gca,'Orientation (degrees)')
    %     ylabel(gca, 'Time (min)')



    reg = regexp(name, 'MMStack_Default.ome', 'split');
    reg = reg{1};
    title(t, strrep(reg,'_', ' ' ));

    saveas(gcf, [outname '_Summary_Plots.png'])



    %% Copy to server
    if uploadresults == 1
        if isremote == 0  % if working with local files, upload to serverfolder (specified in settings)
            [folder2Copy, ~]= fileparts(path);
            lastfolder = regexp(folder2Copy, '\', 'split');
            serverLocation = [serverfolder '\' lastfolder{end} '\'];
            [status,message,messageId]=copyfile(folder2Copy, serverLocation);


        elseif isremote == 1  % if working with remote files, moved analyzed results back to where we found them.
            clear('img')

            tifFiles = dir([tempfolder '\*.tif']);
            for j = 1:length(tifFiles)
                delete(fullfile(tifFiles(j).folder,tifFiles(j).name))
            end

            otherFiles = dir(tempfolder);
            otherFiles = otherFiles(3:end);


            for ri = 1:length(otherFiles)   % copy results to server and clean up our mess.
                file2copy = fullfile(otherFiles(ri).folder, otherFiles(ri).name);
                [uploadLocation, ~]= fileparts(remotepath);
                [status,message,messageId]= copyfile(file2copy,[uploadLocation '\']);
                if status == 1
                    delete(fullfile(otherFiles(ri).folder, otherFiles(ri).name));
                end
            end

        end
    end


    if exist('wormdata', 'var')
        clear('wormdata');
    end

    if exist('img', 'var')
        clear('img')
    end

    if exist('TIFFStack.m','file')
        warning(war);
    end

    %     close all

end