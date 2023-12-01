fld = 'C:\Users\Jeremy\Desktop\230331_zfex813_wt_2'; % Folder containing the data you want to analyze
serverfolder = 'Y:\Calcium Imaging\Wookyu\MuscleCalcium_B12\230331_zfex813_wt_B_2';  % upload everything to this location.

%% settings
startIndex = 1; % which video to start analysis.
startframe =1; % what frame to begin analysis

uploadresults = 0; % upload data to remote location (serverfolder)?
isremote = 0;    % is our tiff file on the server? If so, we'll copy to local
% folder to run the analysis then move the results to the
% server.

plotstuff = 1; % display tracking
videostuff =1; % record video 
framerate = 7; % display video/plots every Nth iteration of loop.
fps = 15;      % frames per sec of input tiff.
troubleshoot =0; % show binary images instead of regular plots
showNormals = 1;
showAxialSignal = 1;

crop = 3; % num pixels to crop from image edge. set to 0 for no cropping.

useautothreshold = 1;% set to 1 to calculate a threshold for each image.
useadaptivethreshold = 0; % if useautothreshold is set to 0 adaptive thresholding can be used
removevignette = 30; % if not zero, size of kernel to use for flatfield correction.

loadtiff =1; % read entire tiff into memory? faster analysis but requires more ram.
minwormarea = 10000; %lower limit to worm area
maxwormarea = 20000; % upper limit to worm area
axSigLen = 200; % how many pixels to use for registering axial signal.
axSigHeight = 20; % how many pixels to use sample perpindicular to the midline.
saveAxialMatrix = 0;

%%
tdir = dir([fld '\**\*.tif']);

for nf =startIndex:length(tdir)
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
        figure('Position',[380.2000 77 690.4000 688.8000],'Color',[1 1 1]);
            tiledlayout(9,3,'TileSpacing', 'compact', 'Padding','tight')
            ax1 = nexttile([3 1]);
            ax2 = nexttile([3 1]);
            ax3 = nexttile([3 1]);
            ax4 = nexttile([2 3]);
            ax7 = nexttile([2 3]);
            areaAx = nexttile([2 3]);
       end


        if videostuff == 1
            vidfig = gcf;
            if exist('v','var') == 1
                close(v)
            end
            videopath = strrep(path, '_MMStack_Default.ome.tif', '_Tracking_Video.mp4');
            v = VideoWriter(videopath,'MPEG-4');
            v.FrameRate = 60;
            open(v)
        end
    end

    %     if showNormals ==1
    %     normfig = figure();
    %     normAx = axes(Parent=normfig);
    %     end

    tic
    info = imfinfo(path);
    toc
    imgWidth = info(1).Width;
    imgHeight = info(1).Height;



    axialSignal = NaN(length(info)/2, axSigLen);
    topAxialSignal = NaN(length(info)/2,axSigLen);
    bottomAxialSignal = NaN(length(info)/2,axSigLen);
    sumSignal = NaN(length(info)/2,1);
    bulkSignal = NaN(length(info)/2,1);
    bulkAboveBkg = NaN(length(info)/2,1);
    backgroundSignal = NaN(length(info)/2,1);
    orientation = NaN(length(info)/2,1);
    area = NaN(length(info)/2,1);
    mag = NaN(length(info)/2,1);
    quadTop = NaN(length(info)/2,1);
    quadBot = NaN(length(info)/2,1);
    time = linspace(0,round((length(info)/2)/fps/60,1),ceil(length(info)/2)); %minutes per frame
    wormIdx = [];

    if saveAxialMatrix == 1
        axialMatrix = NaN(axSigHeight, axSigLen,length(info)/2);
    end


    %% Get comments if available
    %     commentPath = [tdir(nf).folder '\comments.txt'];
    %     if exist(commentPath, 'file')
    %         s = readlines(commentPath);
    %         comments = s(12,:);
    %         tem = regexpi(comments, 'temp', 'split');
    %
    %
    %
    %     end


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
        BW = bwmorph(BW,'fill');


        tempb = BW;
        BW = imdilate(BW,strel('disk',7));
        BW = imerode(BW,strel('disk',7));

        %         BW = imerode(BW,strel('disk',7));
        %         B =bwmorph(B,'open',1);


        tempb2 = BW;

        % identify connected objects
        CC = bwconncomp(BW);
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
                Clen = axSigHeight; % length of perpendicular line to sample

                temptrace = cell(1,length(sortSkel)-1); %NaN(Clen,length(sortSkel)-1);
                tempbf = cell(1,length(sortSkel)-1); %NaN(Clen,length(sortSkel)-1);
                perpX = cell(1,length(sortSkel)-1);
                perpY = cell(1,length(sortSkel)-1);



                parfor ii = 1:length(sortSkel)-1
                    if ii+stepSize<length(sortSkel)
                        seg = sortSkel(ii:ii+stepSize,:);
                    else
                        seg = sortSkel(ii:end,:);
                    end



                    A = [seg(1,2) seg(1,1)]; % [x y] coords of point 1
                    B = [seg(end,2) seg(end,1)]; % [x y] coords of point 2

                    AB = B - A;     % Call AB the vector that points in the direction from A to B

                    % Normalize AB to have unit length
                    AB = AB/norm(AB);

                    % compute the perpendicular vector to the line
                    % because AB had unit norm, so will ABperp
                    ABperp = AB*[0 -1;1 0];

                    % midpoint between A and B
                    ABmid = (A + B)/2;
                    % Compute new points C and D, each at a ditance
                    % Clen off the line. Note that since ABperp is
                    % a vector with unit eEuclidean norm, if I
                    % multiply it by Clen, then it has length Clen.
                    C = ABmid + Clen*ABperp;
                    D = ABmid - Clen*ABperp;

                    try
                        temptrace(ii) = {improfile(GFP,[C(1);D(1)],[C(2);D(2)],Clen)};
                        tempbf(ii) = {improfile(mCh,[C(1);D(1)],[C(2);D(2)],Clen)};
                        perpX(ii) = {[C(1); D(1)]}
                        perpY(ii) = {[C(2); D(2)]}
                    catch
                    end
                end

                hold off



                temptrace = cell2mat(temptrace);
                tempbf = cell2mat(tempbf);
                perpX = cell2mat(perpX);
                perpY = cell2mat(perpY);

                %                 h1 = plot(perpX,perpY);




                if ~isempty(temptrace)
                    % resample axial images
                    x1 = linspace(1,size(temptrace,2), size(temptrace,2));
                    x2 = linspace(1,size(temptrace,2),axSigLen);
                    temptrace = interp1(x1, temptrace',x2)';
                    tempbf = interp1(x1, tempbf',x2)';
                    
                    

                    tt = max(temptrace);


                    % % % % real-time autoFixSignal % % %
                    querryLength = length(tt)*0.1; % fraction of signal to querry
                    leftMean = mean(tt(1:querryLength),'omitnan');
                    rightMean = mean(tt(length(tt)-querryLength:length(tt)),'omitnan');

                    if leftMean<rightMean
                        tt = fliplr(tt);
                        temptrace = rot90(temptrace,2);

                    end

                    top = temptrace(1:size(temptrace,1)/2,:);
                    bot = temptrace(size(temptrace,1)/2+1:end,:);

                     % % % % % % % % % % % % % % % % % % % %


                    axialSignal(i,1:size(tt,2)) = tt;
                    topAxialSignal(i,1:size(top,2)) = max(top);
                    bottomAxialSignal(i,1:size(bot,2)) = max(bot);

                    seg = 57:80;
                    quadTop(i) = mean(max(top(:,seg)));
                    quadBot(i) = mean(max(bot(:,seg)));
                    

                    if saveAxialMatrix == 1
                        axialMatrix(:,:,i) = temptrace;
                    end
                end
            end

            % Bulk signal and background signal
            blksig = GFP(mask);
            sumSignal(i,1) = sum(blksig,"all",'omitnan');
            bulkSignal(i,1) = mean(blksig,"all",'omitnan');
            backgroundSignal(i,1) = mean(GFP(~mask),'all','omitnan');

            abovebkg = blksig>mean(GFP(~mask));
            bulkAboveBkg(i,1) = mean(blksig(abovebkg)-mean(GFP(~mask)),"all",'omitnan');

            area(i,1) = bwprops(wormIdx).Area;



            if plotstuff == 1
                if mod(i,framerate) == 0

                    imshow(label2rgb(L,'jet','k','shuffle'),'Parent', ax1)
                    title(ax1,'Binary Mask');
                    if showNormals == 1
                        line(perpX,perpY,'Color', [0.9 0.9 0.9],'Parent', ax1)
                        title(ax1,'Binary Mask + Normal Vectors');
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
                                gpaTrace = temptrace(:,abs(gdiff):size(GFP, 2));
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
                        axsig = smoothdata(topAxialSignal(1:i,:),1,'gaussian',15)';                        
                        imagesc(axsig,'Parent',ax4)  
                        rectangle(ax4, Position=[1 seg(1) size(topAxialSignal,1), seg(end)-seg(1)], EdgeColor=[0 .45 .74])
                        ax4.CLim = [3000 20000];
                        ax4.XLim = [1, length(topAxialSignal)];
                        ax4.XAxis.Visible = 0;
                        yticks(ax4,[30 size(topAxialSignal,2)-20])
                        yticklabels(ax4,{'Head', 'Tail'})
                        title(ax4,'"Top" muscles')
                        ax4.YAxis.TickLabelRotation = 90;
                        ax4.XAxis.Visible = 0;
                        box(ax4, 'off')
                        ax4.TickLength = [0.005 0.005];
                        colormap turbo

                        
                        
                    end

                    % Bottom signal
                        axsig = smoothdata(bottomAxialSignal(1:i,:),1,'gaussian',15)';
                        imagesc(axsig,'Parent',ax7)  
                        rectangle(ax7,Position=[1 seg(1) size(bottomAxialSignal,1), seg(end)-seg(1)], EdgeColor=[.85 .33 .1])
                        ax7.CLim = [3000 20000];
                        ax7.XLim = [1, length(bottomAxialSignal)];
                        yticks(ax7,[30 size(bottomAxialSignal,2)-20])
                        yticklabels(ax7,{'Head', 'Tail'})
                        title(ax7,'"Bottom" muscles')
                        ax7.YAxis.TickLabelRotation = 90;
                        ax7.XAxis.Visible = 0;
                        box(ax7, 'off')
                        ax7.TickLength = [0.005 0.005];
                        colormap turbo


                    % Area
                    plot(time(1:i),smoothdata(quadTop(1:i), 'gaussian', 3),...
                        time(1:i),smoothdata(quadBot(1:i), 'gaussian', 3),...
                        'Parent', areaAx)
                    ylabel(areaAx, 'quadrent signal');
                    areaAx.TickLength = [0.005 0.005];
                    box(areaAx, 'off');
                    xlim(areaAx,[0 time(end)]);
                    xlabel(areaAx,'Time (min)');

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
    %% %%%%%%%%%%%%%%%%%%%%% Auto Fix Axial Signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fractionToQuerry = 0.1;
    autoAxialSignal = autoFixSignal(axialSignal,fractionToQuerry);

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
    [pk,loc,w] = findpeaks(bulkSignal,'MinPeakProminence',500,'MinPeakDistance',150);
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
    wormdata.topAxialSignal = topAxialSignal;
    wormdata.bottomAxialSignal = bottomAxialSignal;
    if saveAxialMatrix == 1
        wormdata.axialMatrix = axialMatrix;
    end
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
    imagesc(smoothdata(topAxialSignal,1,'gaussian',15)',[3000, 20000]);
    rectangle(Position=[1 seg(1) size(topAxialSignal,1), seg(end)-seg(1)], EdgeColor=[0 .45 .74])
    yticks([30 size(topAxialSignal,2)-20])
    yticklabels({'Head', 'Tail'})
    ylabel('"Top" muscles')
    ax = gca;
    ax.YAxis.TickLabelRotation = 90;
    ax.XAxis.Visible = 0;
    box off



    nexttile([1 1])
    plot(pktime, pktraces, 'Color', [0.7 0.7 0.7])
    hold on
    plot(pktime, pkmean, 'Color', [1 0 0], 'LineWidth', 2);
    hold off
    title(gca, 'Spike Profile');
    ylabel(gca,'Fluorescence (a.u.)');
    colormap bone


    nexttile([1 3]);
    imagesc(smoothdata(bottomAxialSignal,1,'gaussian',15)', [3000, 20000]);
    rectangle(Position=[1 seg(1) size(bottomAxialSignal,1), seg(end)-seg(1)], EdgeColor=[.85 .33 .1])
    yticks([30 size(bottomAxialSignal,2)-20])
    yticklabels({'Head', 'Tail'})
    ylabel('"Bottom" muscles')
    ax = gca;
    ax.YAxis.TickLabelRotation = 90;
    ax.XAxis.Visible = 0;
    box off
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
    ax = gca;
    plot(time,smoothdata(quadTop, 'gaussian', 3),...
        time,smoothdata(quadBot, 'gaussian', 3),...
        'Parent', ax)
    ylabel(ax, 'quadrent signal');
    ax.TickLength = [0.005 0.005];
    box(ax, 'off');
    xlim(ax,[0 time(end)]);
    xlabel(ax,'Time (min)');

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
            clear('img')
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


