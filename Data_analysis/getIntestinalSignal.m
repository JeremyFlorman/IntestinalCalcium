function [bulkSignal,rawAxialSignal,outline,mask,skel,area] = getIntestinalSignal(mCh,gfp,showwork)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
filepath = 'C:\Users\Jeremy\Desktop\220511_zfis178_ser-2;tyra-2;tyra-3-30mM-TA_1\220511_zfis178_ser-2;tyra-2;tyra-3-30mM-TA_1_MMStack_Default.ome.tif';

info = imfinfo(filepath);

downsampleFactor = 1;               % how many lines should we use to sample axial signal. # possible lines/downsample factor.


axialTrace = nan(ceil(length(info)/2), ceil((0.85*info(1).Width)/downsampleFactor));
bulkSignal = nan(ceil(length(info)/2),1);

goparallel = 1;
showwork = 0;
breakfail = 0;

if showwork == 1
    figure();
    ax1 = axes;
    goparallel = 0;
end



%%
tic();


for i = 1:ceil(length(info)/2)
    mCh = imread(filepath,i*2-1);
    gfp = imread(filepath, i*2);

    worm = segWorm(mCh,1,0,0);
    if ~isempty(worm)
        outline = worm.contour.pixels;
        skel = worm.skeleton.pixels;

        bwOutline = zeros(info(1).Height,info(1).Width);
        for k=1:length(outline)
            bwOutline(outline(k,1),outline(k,2)) = 1;
        end

        bwOutline = bwmorph(bwOutline, 'spur',inf);
        bwMask = bwfill(bwOutline,'holes');
        bulkSignal(i) = mean(gfp(bwMask),'all');

        bwSkel = zeros(info(1).Height,info(1).Width);
        for k = 1:length(skel)
            bwSkel(skel(k,1),skel(k,2)) = 1;
        end

        
        se = strel('diamond',2);
        thickskel = imdilate(bwSkel,se);



        bwOutline(logical(thickskel)) = 0;



        %         li1 = ismember(outline,skel(1:10,:),'rows');            % find head/tail (intersection of skeleton and outline)
        %         [rowIndex1,~] = find(li1);
        %         ep1 = outline(rowIndex1, :);                            % coordinates of endpoint 1
        %
        %         li2 = ismember(outline, skel(end-10:end,:),'rows');     % find tail/head (intersection of skeleton and outline)
        %         [rowIndex2,~] = find(li2);
        %         ep2 = outline(rowIndex2, :);                            % coordinates of endpoint 2
        %
        %         bwOutline(ep1(:,1),ep1(:,2)) = 0;  % Split the outline at the endpoints.
        %         bwOutline(ep2(:,1),ep2(:,2)) = 0;


        bwOutline = bwareafilt(bwOutline, 2);


        CC = bwconncomp(bwOutline,8);
        if CC.NumObjects == 1
            breakfail = breakfail+1;
        end

        if CC.NumObjects >1

            [flank1R, flank1C] = ind2sub(size(bwOutline), CC.PixelIdxList{1,1});
            [flank2R, flank2C] = ind2sub(size(bwOutline), CC.PixelIdxList{1,2});


            if length(flank1R)>20 && length(flank2R)>20

                P1 = neighborSearch(bwOutline,flank1R,flank1C);
                P2 = neighborSearch(bwOutline,flank2R,flank2C);

                if ~isempty(P1) && ~isempty(P2)

                    trace1 = Flail(bwOutline,P1(1,:),inf,'clockwise');
                    trace2 = Flail(bwOutline,P2(1,:),inf,'clockwise');

                    %%

                    if showwork == 1
                        imshow(gfp,[0 20000],'Parent', ax1)
                        hold on

                    end

                    shortertrace = min(length(trace1),length(trace2));
                    if goparallel == 0
                        profile = nan(ceil(shortertrace/(2*downsampleFactor)),10);

                        for k = 1:ceil(shortertrace/(2*downsampleFactor))
                            profile(k,1:10) = improfile(gfp,[trace1(k,2) trace2(k,2)],[trace1(k,1), trace2(k,1)],10);

                            if showwork == 1
                                plot([trace1(k,2) trace2(k,2)],[trace1(k,1), trace2(k,1)],'Parent', ax1);
                            end
                        end

                    elseif goparallel == 1
                        profile = nan(ceil(shortertrace/(2*downsampleFactor)),1);
                        parfor k = 1:ceil(shortertrace/(2*downsampleFactor))
                            profile(k) = max(improfile(gfp,[trace1(k,2) trace2(k,2)],[trace1(k,1), trace2(k,1)],10));
                        end
                    end

                    if goparallel == 0
                        mp = max(profile,[],2);
                    elseif goparallel == 1
                        mp = profile;
                    end

                    mp = resample(mp, size(axialTrace,2),length(mp),5,20);
                    axialTrace(i,1:length(mp)) = mp;
                    %         plot(mp)
                    if showwork == 1
                        drawnow
                        hold off
                    end
                end
            end
        end
    end
    if mod(i,90) == 0
        imagesc(axialTrace,[0 50000])
        colormap turbo
        drawnow
    end

    if mod(i,90) == 0
        disp([num2str((i/ceil(length(info)/2))*100) ' % done in ' num2str(toc()) 'seconds'])
    end
end

axialTrace = autoFixSignal(axialTrace);
figure('Position', [661.8000 65 378.4000 673.600], 'Color',[1 1 1])
imagesc(axialTrace, [0 50000])
%%
wormdata.autoAxialSignal = axialTrace;
wormdata.bulkSignal = bulkSignal;


datasavename = strrep(filepath, 'MMStack_Default.ome.tif', 'wormdata.mat');
save(datasavename, 'wormdata')


end