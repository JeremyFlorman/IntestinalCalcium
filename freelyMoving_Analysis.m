

fld = 'E:\Jeremy Acquisitions\DMP_Mutants\gon-2(kd)' % Folder containing the data you want to analyze
serverfolder = 'Z:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\gon-2(kd)'; % upload everything to this location.

uploadresults = 1; % upload data to remote location (serverfolder)?
isremote = 0; % is our tiff file on the server? If so, we'll copy to local
% folder to run the analysis then move the results to the
% server.
 

%% settings
plotstuff = 1; % display tracking
videostuff =1; % record video
framerate = 15; % display video/plots every Nth iteration of loop.
troubleshoot =0; % show binary images instead of regular plots
showtranslation = 0;

crop = 7; % num pixels to crop from image edge. set to 0 for no cropping.
minwormarea = 10000; %lower limit to worm area
maxwormarea = 20000; % upper limit to worm area
useautothreshold = 1;% set to 1 to calculate a threshold for each image.
useadaptivethreshold = 0; % if useautothreshold is set to 0 adaptive thresholding can be used
removevignette = 160; % if not zero, size of kernel to use for flatfield correction
bearing = 1;   % sets the initial direction of the head (1 or 0).
startframe =1; % when to begin analysis
loadtiff =1; % read entire tiff into memory? faster analysis but requires more ram.

%%
tdir = dir([fld '\**\*.tif']);

for nf = 1:length(tdir)
    path = fullfile(tdir(nf).folder, tdir(nf).name)


    if isremote == 1
        tic
        tempfolder = 'E:\tmp'; %set temp directory for copying tiff files.
        copyfile(path, tempfolder)
        remotepath = path;                         % save the original path for later uploading.
        path = fullfile(tempfolder, tdir(nf).name); % rename path to the local path.
        toc
    end

    %%
    if plotstuff == 1
        figure('Position', [506 92 1332 849],'Color',[1 1 1]);
        tiledlayout(4,3)
        ax1 = nexttile([2 1]);
        ax2 = nexttile([2 1]);
        ax3 = nexttile([2 1]);
        ax4 = nexttile([1 1]);
        ax5 = nexttile([1 1]);
        ax6 = nexttile([1 1]);
        ax7 = nexttile([1 3]);


        if videostuff == 1
            vidfig = gcf;
            if exist('v','var') == 1
                close(v)
            end
            videopath = strrep(path, '_MMStack_Default.ome.tif', '_Tracking_Video.mp4');
            v = VideoWriter(videopath,'MPEG-4');
            v.FrameRate = 5;
            open(v)
        end
    end

    if showtranslation ==1
        translationfig = figure();
        transax = axes('Parent',translationfig);
    end

    tic
    info = imfinfo(path);
    toc
    imgWidth = info(1).Width;
    imgHeight = info(1).Height;

    axialSignal = NaN(length(info)/2, max([imgWidth imgHeight]));
    autoAxialSignal = NaN(length(info)/2, max([imgWidth imgHeight]));
%     autoAxialBF = NaN(length(info)/2, max([imgWidth imgHeight]));
    axialBF = NaN(length(info)/2, max([imgWidth imgHeight]));
%     axmat = cell(length(info)/2,2);
    sumSignal = NaN(length(info)/2,1);
    bulkSignal = NaN(length(info)/2,1);
    bulkAboveBkg = NaN(length(info)/2,1);
    backgroundSignal = NaN(length(info)/2,1);
    orientation = NaN(length(info)/2,1);
    area = NaN(length(info)/2,1);
    mag = NaN(length(info)/2,1);
    fps = 15;
    time = linspace(0,ceil((length(info)/2)/fps/60),ceil(length(info)/2)); %minutes per frame
    wormIdx = [];

    if loadtiff == 1
        if exist('TIFFStack.m','file')
        img = TIFFStack(path);
         war = warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
         warning('off', 'MATLAB:imagesci:tifftagsread:expectedTagDataFormat');
         warning('off','imageio:tiffmexutils:libtiffWarning')
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
        %imshow(imadjust(mCh))
        % set up thresholding and binary image processing
        if useautothreshold ==1
            T = graythresh(mCh);
            B = imbinarize(mCh, T); % create mask
        elseif useadaptivethreshold == 1
            B = imbinarize(mCh,'adaptive','ForegroundPolarity','dark','Sensitivity',0.3);
        end



        se = strel('disk',7);
        B = imcomplement(B);

        tempb = B;
        B = imdilate(B,se);
        B = imerode(B,se);

        tempb2 = B;

        % identify connected objects
        CC = bwconncomp(B);
        L = labelmatrix(CC);
        bwprops = regionprops(L, 'Area', 'Centroid','Orientation');

        % Filter connected components to get biggest most central object.
        xy = vertcat(bwprops.Centroid);
        x = xy(:,1);
        y = xy(:,2);
        distances = sqrt((imgWidth/2 - x) .^ 2 + (imgHeight/2 - y) .^ 2);
        [centralSize, centralIdx] = min(distances); % most central object
        [bigSize, bigIdx] = max([bwprops.Area]); % largest object

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
            wormIdx = [];
        end

        % create a copy of the label matrix Lw that contains only the worm.
        Lw = L;
        %   imshow(label2rgb(L,'jet','k','shuffle'))
        if ~isempty(wormIdx)
            Lw(Lw~=wormIdx) = 0;



            % generate mask, outline and skeleton
            mask = logical(Lw);
            outline = bwmorph(mask, 'remove',1);
            %   skel = bwskel(mask,'MinBranchLength', 15);
            skel = bwmorph(mask,'thin', inf);
            outskel = logical(outline+skel);


            % % Oreintation Block: orientation is determined from the regoion props
            % % orientation value of the mask. this value is stored in the array
            % % 'orientation'. This stored record will be used to detect flips in
            % % oreintation (from -90 to 90 or vice versa). A second variable
            % % 'rotation' is used to rotate images to keep worms horizontal. By
            % % modifying the sign of 'rotation' based on detection of flipping events
            % % it should be possible to keep the worm facing the same direction over
            % % the course of the experiment.The variable bearing is defiend in the
            % % settings and can be used to flip the initial orientation of the worm.

            orientation(i,1) = bwprops(wormIdx).Orientation;

            if ~isnan(orientation(i,1))
                rotation = 180-orientation(i,1);
            else
                rotation = 0;
            end


            % rotate images relative to major axis of worm.
            rCherry = imrotate(mCh, rotation,'bicubic', 'loose');
            rGFP = imrotate(GFP, rotation,'bicubic', 'loose');
            rOutskel = imrotate(outskel, rotation, 'loose');
            rMask = imrotate(mask, rotation, 'loose');
            
            rSkel = bwskel(rMask,'MinBranchLength', 15);
            %             rSkel = bwmorph(rMask,'thin', inf);
            [ep] = bwmorph(rSkel,'endpoints');

            if nnz(ep) >0
                [ey,ex] = find(ep,1);
                sortSkel= bwtraceboundary(rSkel,[ey ex],'E');
                sortSkel = sortSkel(1:ceil(length(sortSkel)/2),:);
                %             nnz(ep)
                %Scanning Block: since the worm midline defined by the skeleton does
                %not always hit the portion of the intestine we want to measure it is
                %useful to scan multiple times. This block translates the midline over
                %a 15 pixel distance (as defined by ii) across the y axis. This relies
                %on the fact that the worm is horizontal due to rotation. The values of
                %each scanning iteration are collected in 'temptrace' and the median is
                %resampled into 'axialSignal'

                temptrace = NaN(11,length(sortSkel));
                tempbf = NaN(11,length(sortSkel));

                for ii = -10:10

                    skelx = sortSkel(:,2);
                    skely = sortSkel(:,1)+ii;




                    if showtranslation == true
                        % old tracing
                        tempskel = imtranslate(rSkel, [0, ii]);
                        [oldy, oldx] = find(tempskel);
                        if nnz(tempskel)<= length(temptrace)
                            temptrace(ii+11,1:nnz(tempskel)) = rGFP(tempskel);
                            tempbf(ii+11,1:nnz(tempskel)) = rCherry(tempskel);
                        end

                        imshow(imoverlay(imadjust(rCherry), tempskel, [1 0 0]), 'Parent', transax)
                        hold on
                        for oldidx = 1:length(oldy)
                            scatter(oldx(oldidx),oldy(oldidx),30, [0.1 0.1 1])
                            drawnow()
                        end
                        hold off
                        % new tracing
                        imshow(imoverlay(imadjust(rCherry), ep, [1 0 0]), 'Parent', transax)
                        hold on
                        for tidx = 1:length(skely)
                            scatter(skelx(tidx),skely(tidx),30)
                            drawnow()
                        end
                        hold off


                    end

                    try
                        for fidx = 1:length(skelx)
                            temptrace(ii+11, fidx) = rGFP(skely(fidx), skelx(fidx));
                            tempbf(ii+11, fidx) = rCherry(skely(fidx), skelx(fidx));
                        end
                    catch
                    end

                end

                if ~isempty(temptrace)
                    tt = resample(max(temptrace), size(axialSignal,2), size(temptrace,2),5,20);  % max?
                    abf = resample(max(tempbf), size(axialSignal,2), size(tempbf,2),5,20);

%                     axmat(i,1) = {tempbf};
%                     axmat(i,2) = {temptrace};
                    axialBF(i,1:size(abf,2)) = abf;
                    axialSignal(i,1:size(tt,2)) = tt;
                end
            end
            
%% Bulk signal and background signal
            blksig = rGFP(rMask);
            sumSignal(i,1) = sum(blksig,"all",'omitnan');
            bulkSignal(i,1) = mean(blksig,"all",'omitnan');
            backgroundSignal(i,1) = mean(GFP(~mask),'all','omitnan');
            
            abovebkg = blksig>mean(GFP(~mask));
            bulkAboveBkg(i,1) = sum(blksig(abovebkg)-mean(GFP(~mask)));
            


            area(i,1) = bwprops(wormIdx).Area;

%%

            if plotstuff == 1
                if mod(i,framerate) == 0
                    %                     imshow(axmat{i,1}, [1000 12000], 'Parent', ax1)
                    %                     colormap bone
                    imshow(label2rgb(L,'jet','k','shuffle'),'Parent', ax1)
                    title(ax1,'Binary Masks');

                    if troubleshoot == 1
                        imshow(tempb,'Parent', ax2);
                        title(ax2,'Initial Threshold');

                        imshow(tempb2,'Parent', ax3);
                        title(ax3,'Processed Mask');
                    else
                        try
                            mdiff = size(rCherry,2)-size(tempbf,2);
                            mpadTrace = padarray(tempbf,[0, ceil(mdiff/2)],0,'both');
                            mpadTrace = mpadTrace(:,1:size(rCherry, 2));
                            mpad_rOutskel = padarray(rOutskel, [size(mpadTrace,1),0], 'post');
                            mmergedImage = vertcat(rCherry, mpadTrace);
                            mmergedOverlay = imoverlay(imadjust(mmergedImage, [0 0.2]), mpad_rOutskel, [1 0 0]);
                            imshow(mmergedOverlay,'Parent', ax2)
                            title(ax2,'Brightfield');

                            gdiff = size(rGFP,2)-size(temptrace,2);
                            gpadTrace = padarray(temptrace,[0, ceil(gdiff/2)],0,'both');
                            gpadTrace = gpadTrace(:,1:size(rGFP, 2));
                            gpad_rOutskel = padarray(rOutskel, [size(gpadTrace,1),0], 'post');
                            gmergedImage = vertcat(rGFP, gpadTrace);
                            gmergedOverlay = imoverlay(imadjust(gmergedImage, [0.06 0.2]), gpad_rOutskel, [0 1 0]);
                            imshow(gmergedOverlay,'Parent', ax3)
                            title(ax3,'GCaMP');
                        catch
                            imshow(imoverlay(imadjust(rCherry, [0 0.2]), rOutskel, [1 0 0]), 'Parent', ax2)
                            title(ax2,'Brightfield');

                            imshow(imoverlay(imadjust(rGFP,[0.06 0.2]), rOutskel, [0 1 0]), 'Parent', ax3)
                            title(ax3,'GCaMP');
                        end
                    end

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


                    plot(time,bulkSignal(:),time,backgroundSignal(:), 'Parent', ax7)
                    title(ax7, 'Whole Body Ca^2^+ Signal')
                    ylabel(ax7, 'Mean Fluorescent Intensity (a.u.)');
                    xlabel(ax7,'Time (min)');
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

    disp('file processed in:')
    toc

    if exist('v','var') == 1
        close(v)
    end
    %%%%%%%%%%%%%%%%%%%%%%% Auto Fix Axial Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:length(axialSignal)
        left = mean(axialSignal(ii,1:30),'omitnan');
        right = mean(axialSignal(ii,end-30:end),'omitnan');
        if left > right
            autoAxialSignal(ii,:) = fliplr(axialSignal(ii,:));
%             autoAxialBF(ii,:) = fliplr(axialBF(ii,:));
%             axmat(ii,1) = {fliplr(axmat{ii,1})};
%             axmat(ii,2) = {fliplr(axmat{ii,2})};
        elseif left <= right
            autoAxialSignal(ii,:) = axialSignal(ii,:);
%             autoAxialBF(ii,:) = axialBF(ii,:);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    %% peak analysis
    [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',500,'MinPeakDistance',150);
    peakpad = 15*15; % framerate*time in seconds;
    pktime = linspace(-15,15, peakpad*2)';
    if isempty(loc)
        loc = NaN;
        locinmin = NaN;
    else
        locinmin = loc/(15*60);
    end

    pktraces = NaN(peakpad*2,length(loc));
    for k = 1:length(loc)
        prepeak = loc(k)-peakpad;
        postpeak = loc(k)+peakpad-1;
        temppeak = [];

        if prepeak>1 && postpeak<length(bulkSignal)
            pktraces(1:peakpad*2,k) = bulkSignal(prepeak:postpeak,1);
        end

    end
    pkmean = mean(pktraces,2,'omitnan');

    %% Save Stuff
    [folder, name, ext] = fileparts(path);
    outname = fullfile(folder, name);
    datasavename = strrep(outname, 'MMStack_Default.ome', 'wormdata.mat');


    wormdata = struct();
    wormdata.autoAxialSignal = autoAxialSignal;
    wormdata.sumSignal = sumSignal;
    wormdata.bulkSignal = bulkSignal;
    wormdata.bulkAboveBkg = bulkAboveBkg;
    wormdata.backgroundSignal = backgroundSignal;
    wormdata.orientation = orientation;
    wormdata.area = area;
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
        time = linspace(0,ceil((length(info)/2)/fps/60),ceil(length(info)/2)); %minutes per frame
    end
    if ~exist('pk','var')
        [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',1000, 'MinPeakDistance',150);
        peakpad = 15*15;
        pktime = linspace(-15,15, peakpad*2)';
        pkmean = mean(pktraces,2,'omitnan');
    end



    figure('Position', [135.4000 142.6000 1232 586.4000])
    t = tiledlayout(3,4);

    nexttile([1 3])
    if ~isnan(locinmin)
        plot(time,bulkSignal,locinmin,pk*1.01, 'rv')
    else
        plot(time,bulkSignal)
    end

    hold on
    plot(time, backgroundSignal)
    hold off

    xlim([0 10])
    xlabel(gca, 'Time (min)')
    ylabel(gca,'Fluorescence (a.u.)');
    title(gca, 'Whole Animal Calcium Trace')


    ax = nexttile([3 1]);
    imagesc(smoothdata(autoAxialSignal,1,'gaussian',60))
    title(gca, 'Axial Calcium Trace')
    ct = gca;
    hold on
    plot(size(autoAxialSignal,2),loc, '<w', 'MarkerFaceColor' ,[.4 .5 .6]);
    hold off


    nexttile([1 1])
    plot(pktime, pktraces, 'Color', [0.7 0.7 0.7])
    hold on
    plot(pktime, pkmean, 'Color', [1 0 0], 'LineWidth', 2);
    hold off
    title(gca, 'Spike Profile');
    ylabel(gca,'Fluorescence (a.u.)');
    colormap bone

    nexttile([1 1])
    histogram(w./15,5);
    ylim([0 10])
    xlim([0 15])
    title(gca,'Peak Widths');
    ylabel(gca,'Count');
    xlabel(gca,'Time (s)');

    nexttile([1 1])
    histogram(diff(loc)./15, 5);
    title(gca,'Inter-Peak Interval');
    ylim([0 10])
    xlim([0 120])
    xlabel(gca,'Time (s)');
    ylabel(gca,'Count');

    nexttile([1 1]);
    plot(orientation,time);
    ylim([0 10])
    title(gca, 'Worm Orientation')
    set(gca,'YDir', 'reverse')
    xlabel(gca,'Orientation (degrees)')
    ylabel(gca, 'Time (min)')


    nexttile([1 1]);
    plot(time,area)
    xlim([0 10])
    title(gca, 'Worm Area')
    ylabel(gca,'Pixels')
    xlabel(gca,'Time (min)')

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
            pattern = fullfile(tempfolder, '*.tif'); % identify large tif file which was moved to local drive
            tifffile = dir(pattern);
            ld = dir(tempfolder);
            ld = ld(3:end);

            tifflist = zeros(length(ld),1);
            for i = 1:length(ld)
                tifflist(i) =  strcmpi(ld(i).name, tifffile.name);
            end

            noTiff = ld(~tifflist);     % exclude tif file, we dont need to move it to the server.

            for ri = 1:length(noTiff)   % copy results to server and clean up our mess.
                file2copy = fullfile(noTiff(ri).folder, noTiff(ri).name);
                [uploadLocation, ~]= fileparts(remotepath);
                [status,message,messageId]= copyfile(file2copy,[uploadLocation '\']);
                if status == 1
                    delete(fullfile(noTiff(ri).folder, noTiff(ri).name));
                end
            end
            delete(fullfile(tifffile.folder, tifffile.name));
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











%% Move everything to z drive
% localfolders = dir('E:\Jeremy Acquisitions\FreelyMoving\Analysis\');
%
% for i = 3:length(localfolders)
%     query = ['Z:\Calcium Imaging\FreelyMoving\**\*' localfolders(i).name];
%     remotefolder = dir(query);
%     remotepath = fullfile(remotefolder.folder, remotefolder.name);
%
%     localdir = dir(fullfile(localfolders(i).folder, localfolders(i).name));     % move everything
%     ftm = fullfile({localdir.folder},{localdir.name});  %move everything
%     ftm = ftm(4:end);
%     for j = 1:length(ftm)
%          copyfile(ftm{j}, [remotepath '\'])
%      disp(['Moving: ' ftm{j} ' to ' remotepath])
%     end
%
%     axial = dir([fullfile(localfolders(i).folder, localfolders(i).name) '\*axialSignal.txt']); % move axial
%     plots = dir([fullfile(localfolders(i).folder, localfolders(i).name) '\*Summary_Plots.png']);  % move summary plot
%     copyfile(fullfile(axial.folder, axial.name), [remotepath '\'])
%     copyfile(fullfile(plots.folder, plots.name), [remotepath '\'])
% end


