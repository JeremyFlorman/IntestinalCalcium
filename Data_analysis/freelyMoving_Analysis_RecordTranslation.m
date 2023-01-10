% % % Path Selection
% [fn,fp] = uigetfile('E:\Jeremy Acquisitions\FreelyMoving\Analysis\*.tif');
% path = fullfile(fp,fn);
% %
%

%%
%fld = 'E:\Jeremy Acquisitions\Exogenous_TA\5_min_acclimation\'; % do as loop?
fld = 'C:\Users\Jeremy\Desktop\220316_zfis178_wildtype_1' % Folder containing the data you want to analyze
serverfolder = 'Z:\Calcium Imaging\Exogenous_Tyramine\Receptor_Mutants'; % upload everything to this location.

isremote = 0; % is our tiff file on the server? If so, we'll copy to local 
              % folder to run the analysis then move the results to the
              % server. 


%% settings
plotstuff = 1; % display tracking
videostuff =0; % record video
framerate = 15; % display video/plots every Nth iteration of loop.
troubleshoot =0; % show binary images instead of regular plots
showtranslation = 1;
if showtranslation ==1
    if exist('tv','var') == 1
        close(tv)
    end
    transvideopath = [fld '\Translation_Video.mp4'];
    tv = VideoWriter(transvideopath, 'MPEG-4');
    tv.FrameRate = 170;
    open(tv)
end

crop = 7; % num pixels to crop from image edge. set to 0 for no cropping.
minwormarea = 1000; %lower limit to worm area (1000 for adult w/ 5x objective)
maxwormarea = 5000; % upper limit to worm area (5000 for adult w/ 5x objective)
useautothreshold = 1;% set to 1 to calculate a threshold for each image.
useadaptivethreshold = 0; % if useautothreshold is set to 0 adaptive thresholding can be used
removevignette = 80; % if not zero, size of kernel to use for flatfield correction
bearing = 1;   % sets the initial direction of the head (1 or 0).
startframe =1; % when to begin analysis


%%
tdir = dir([fld '\**\*.tif']);

for nf = 1:length(tdir)
    path = fullfile(tdir(nf).folder, tdir(nf).name)

tic
    if isremote == 1
        tic
        tempfolder = 'C:\tmp'; %set temp directory for copying tiff files.
        copyfile(path, tempfolder, 'f')
        remotepath = path;                         % save the original path for later uploading.
        path = fullfile(tempfolder, tdir(nf).name); % rename path to the local path.
    end
    toc
    %%
    if plotstuff == 1
        figure('Position', [506 92 1332 849]);
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
            videopath = strrep(path, '_MMStack_Default.ome.tif', '_Tracking_Video.avi');
            v = VideoWriter(videopath);
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
    axialBF = NaN(length(info)/2, max([imgWidth imgHeight]));
    bulkSignal = NaN(length(info)/2,1);
    orientation = NaN(length(info)/2,1);
    rotation = NaN(length(info)/2,1);
    area = NaN(length(info)/2,1);
    mag = NaN(length(info)/2,1);
    time = linspace(0,10,9000);
    wormIdx = [];
    isflipped = 1;
    %% Tracking Block

    for i = 157:7:857%startframe:length(info)/2
        mCh = imread(path, (i*2)-1);      % mCherry signal is odd # frames?
        GFP = imread(path, i*2);    % GCaMP signal is even # frames?

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
            B = imbinarize(mCh,'adaptive','ForegroundPolarity','dark','Sensitivity',0.5);
        end



        se = strel('disk',3);
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

            if i==startframe
                rotation(i,1)=orientation(i,1);
            elseif i>1
                prevRot = rotation(i-1,1); % find the rotation from the previous frame.

                if isnan(prevRot)  % if the last frame was dropped, search for the most recent previous orientation.
                    pr = find(~isnan(rotation(:,1)),1,'last');
                    if ~isnan(pr)
                        prevRot = rotation(pr,1);
                    else
                        prevRot = 0;
                    end
                end

                if prevRot >0 && orientation(i,1) <=90 % positive rotation no flip
                    rotation(i,1)=orientation(i,1);
                elseif prevRot <0 && orientation(i,1) >=-90 % negative rotation no flip
                    rotation(i,1)=orientation(i,1);
                elseif prevRot == 0
                    rotation(i,1)=0;
                end
            end

            if bearing == 1
                rotation(i,1) = rotation(i,1)*-1+180;
            elseif bearing == 0
                rotation(i,1) = rotation(i,1)*-1;
            end


            % rotate images relative to major axis of worm.
            rCherry = imrotate(mCh, rotation(i,1), 'loose');
            rGFP = imrotate(GFP, rotation(i,1), 'loose');
            rOutskel = imrotate(outskel, rotation(i,1), 'loose');
            %         rSkel = imrotate(skel, rotation(i,1), 'loose');
            rMask = imrotate(mask, rotation(i,1), 'loose');
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
                %a 10 pixel distance (as defined by ii) across the y axis. This relies
                %on the fact that the worm is horizontal due to rotation. The values of
                %each scanning iteration are collected in 'temptrace' and the median is
                %resampled into 'axialSignal'

                temptrace = NaN(11,length(sortSkel));
                tempbf = NaN(11,length(sortSkel));

                for ii = 0 %-5:5

                    skelx = sortSkel(:,2);
                    skely = sortSkel(:,1)+ii;

                    % old tracing protocol
                    tempskel = imtranslate(rSkel, [0, ii]);
                    [oldy, oldx] = find(tempskel);
                    if nnz(tempskel)<= length(temptrace)
                        temptrace(ii+6,1:nnz(tempskel)) = rGFP(tempskel);
                        tempbf(ii+6,1:nnz(tempskel)) = rCherry(tempskel);
                    end


                    if showtranslation == true
                        % old tracing
%                         imshow(imoverlay(imadjust(rCherry), tempskel, [1 0 0]), 'Parent', transax)
%                         hold on
%                         for oldidx = 1:length(oldy)
%                             scatter(oldx(oldidx),oldy(oldidx),30, [0.1 0.1 1])
%                             drawnow()
%                         end
%                         hold off
                        % new tracing
%                         imshow(imoverlay(imadjust(rCherry), ep, [1 0 0]), 'Parent', transax)
                        imshow(imoverlay(imadjust(rGFP), tempskel, [1 0 0]), 'Parent', transax)
                        
                        hold on
                        for tidx = 1:length(skely)
                            scatter(skelx(tidx),skely(tidx),10)
                            
                            drawnow()
                            tf = getframe(translationfig);
                            writeVideo(tv,tf)

                        end
                        


                    end

hold off
                    for fidx = 1:length(skelx)
                        temptrace(ii+6, fidx) = rGFP(skely(fidx), skelx(fidx));
                        tempbf(ii+6, fidx) = rCherry(skely(fidx), skelx(fidx));
                    end

                    
                end

                if ~isempty(temptrace)
                    tt = resample(max(temptrace), size(axialSignal,2), size(temptrace,2));  % max?
                    abf = resample(median(tempbf,1), size(axialSignal,2), size(tempbf,2));


                    axialBF(i,1:size(abf,2)) = abf;
                    axialSignal(i,1:size(tt,2)) = tt;
                end
            end

            if i~= startframe
                prevBF = axialBF(i-1,:);
                currBF = axialBF(i,:);

                samediff = sum(abs(currBF-prevBF));
                flipdiff = sum(abs(fliplr(currBF)-prevBF));

                if flipdiff<samediff
                    isflipped = isflipped*-1;
                end

                if isflipped == -1
                    axialSignal(i,:) = fliplr(axialSignal(i,:));
                end
            end



            bulkSignal(i,1) = mean(mean(rGFP(rMask)));
            area(i,1) = bwprops(wormIdx).Area;

            if plotstuff == 1
                if mod(i,framerate) == 0
                    imshow(label2rgb(L,'jet','k','shuffle'),'Parent', ax1)
                    title(ax1,'Binary Masks');

                    if troubleshoot == 1
                        imshow(tempb,'Parent', ax2);
                        title(ax2,'Initial Threshold');

                        imshow(tempb2,'Parent', ax3);
                        title(ax3,'Processed Mask');
                    else
                        imshow(imoverlay(imadjust(rCherry), rOutskel, [1 0 0]), 'Parent', ax2)
                        title(ax2,'Brightfield');

                        imshow(imoverlay(imadjust(rGFP,[0.05 0.3]), rOutskel, [0 1 0]), 'Parent', ax3)
                        title(ax3,'GCaMP');
                    end

                    plot(time,area(:), 'Parent', ax4);
                    title(ax4,'Worm Area');
                    ylabel(ax4,'Pixels');
                    xlabel(ax4,'Time (min)');

                    plot(time,orientation(:),time,rotation(:), 'Parent', ax5);
                    title(ax5,'Worm Orientation');
                    ylabel(ax5,'Degrees');
                    xlabel(ax5,'Time (min)');

                    plot(1:size(axialSignal,2),axialSignal(i,:),1:size(axialBF,2),axialBF(i,:), 'Parent', ax6)
                    ax6.XLim = [0 size(axialSignal,2)];
                    ax6.YLim = [0 50000];
                    title(ax6,'Signal Along Midline');
                    ylabel(ax6,'Mean Fluorescent Intensity (a.u.)');
                    xlabel(ax6, 'head <--- Distance (pixels) ---> tail');
                    legend(ax6,{'GCaMP6', 'Brightfield'}, 'Location', 'bestoutside');

                    plot(time,bulkSignal(:), 'Parent', ax7)
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
        disp(['Working... ' num2str((i/length(info))*200) '% complete, just chill...'])
    end

close(tv)

    %%%%%%%%%%%%%%%%%%%%%%% Auto Fix Axial Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:length(axialSignal)
        left = mean(axialSignal(ii,1:30),'omitnan');
        right = mean(axialSignal(ii,end-30:end),'omitnan');
        if left > right
            autoAxialSignal(ii,:) = fliplr(axialSignal(ii,:));
        elseif left <= right
            autoAxialSignal(ii,:) = axialSignal(ii,:);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    %% peak analysis
    % figure()
    % findpeaks(bulkSignal,'MinPeakProminence',500, 'MinPeakDistance',150);
    [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',750,'MinPeakDistance',150);
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
            %     elseif prepeak<1
            %         temppeak = bulkSignal(1:postpeak,1);
            %         pktraces(1:length(temppeak),k) = temppeak;
            %     elseif postpeak >length(bulkSignal)
            %         temppeak = bulkSignal(prepeak:end,1);
            %         pktraces(1:length(temppeak),k) = temppeak;
        end

    end
    pkmean = mean(pktraces,2,'omitnan');

    %% Save Stuff
    [folder, name, ext] = fileparts(path);
    outname = fullfile(folder, name);
    datasavename = strrep(outname, 'MMStack_Default.ome', 'wormdata.mat');


    wormdata = struct();
    wormdata.rawAxialSignal = axialSignal;
    wormdata.autoAxialSignal = autoAxialSignal;
    wormdata.axialBrightField = axialBF;
    wormdata.bulkSignal = bulkSignal;
    wormdata.orientation = orientation;
    wormdata.area = area;
    wormdata.peakTraces = pktraces;
    wormdata.peakLoc = loc;
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
        time = linspace(0,10,9000);
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
    xlim([0 10])
    xlabel(gca, 'Time (min)')
    ylabel(gca,'Fluorescence (a.u.)');
    title(gca, 'Whole Animal Calcium Trace')


    ax = nexttile([3 1]);
    imagesc(autoAxialSignal)
    title(gca, 'Axial Calcium Trace')
    ct = gca;
    hold on
    plot(size(axialSignal,2),loc, '<w', 'MarkerFaceColor' ,[.4 .5 .6]);
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


