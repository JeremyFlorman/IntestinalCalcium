folder = 'C:\Users\Jeremy\Desktop\mCherryOP50\220316_zfis178_mcherryOP50_1';
tdir = dir([folder '\*.tif']);
mdir = dir([folder '\*.mat']);

tiff = fullfile(tdir(1).folder, tdir(1).name);
mfile = fullfile(mdir(1).folder, mdir(1).name);

load("rmap2.mat");
load("gmap.mat");

% rmap = rmap.rmap;
% gmap = gmap.gmap;

info = imfinfo(tiff);
%%
w = load(mfile, 'wormdata');
d = w.wormdata;

ypad = 0.01;
scaledata = 1;
recordvid = 1;
timepre = 5;
timepost =6;


for dn = 2%1:length(d.peakLoc(:))

    startidx = d.peakLoc(dn)-(timepre*15);
    endidx = d.peakLoc(dn) + (timepost*15);


savename = [folder '\dmp_' num2str(dn) '.mp4'];




bf = d.fixedAxialBrightField;
gc = d.autoAxialSignal;

if endidx > length(bf)
    endidx = length(bf);
end

mbf = min(bf(startidx:endidx,:),[],"all","omitnan");
xbf = max(bf(startidx:endidx,:),[],"all","omitnan");

mbf = mbf-mbf*ypad;
xbf = xbf+xbf*ypad;

mgc = min(gc(startidx:endidx,:),[],'all',"omitnan");
xgc = max(gc(startidx:endidx,:),[],'all',"omitnan");

mgc = mgc-mgc*ypad;
xgc = xgc+xgc*ypad;

x= 1:size(gc,2);


%%

if exist('v', 'var')
    close(v)
end

vidfig = figure('Position', [504.2000 121 846.4000 652], 'Color', [1 1 1]);

if recordvid == 1
    v = VideoWriter(savename, 'MPEG-4');
    v.Quality = 100;
    v.FrameRate = 15;
    open(v);
end


tiledlayout(2,2, 'TileSpacing','tight')
time = timepre*-1;

t1 = nexttile([1,1]);
t2 = nexttile([1,1]);
t3 = nexttile([1,2]);




for i = startidx:endidx-1
    bfimg = imread(tiff,i*2+1);
    gcimg = imread(tiff,i*2);

    imshow(bfimg,'Colormap', rmap,'DisplayRange', [5000 13000],'Parent',t1)
    time = time+0.0667;
    text(t1,5, 10, [num2str(time) ' sec'], 'FontSize',10)
 
    
    title('mCherry (food)','Parent',t1)
    imshow(gcimg,'Colormap', gmap, 'DisplayRange',[5000 16000],'Parent',t2)
    title('GCaMP (Ca^2^+)','Parent',t2)

    if scaledata == 0
        yyaxis left
        plot(x,bf(i,:),'r','LineWidth',2,'Parent',t3)
        ylim([mbf xbf])

        yyaxis right
        plot(x,gc(i,:),'g','LineWidth',2,'Parent',t3)
        ylim([mgc xgc])

        drawnow()
    elseif scaledata == 1
        bf = rescale(bf);
        gc = rescale(gc);
        plot(x,bf(i,:),'r', x, gc(i,:), 'g', 'LineWidth',2,'Parent',t3);
        ylim([-0.15, 1.15])

    end

    legend({'mCherry(Food)', 'GCaMP(Ca^2^+)'},'Location', 'northwest')
    xlim([-5 261])
    ylabel('Normalized Signal')
    t3.XTick = [20 130 240];
    t3.XTickLabel = {'Head', 'Mid-body', 'Tail'};
%     time = time+0.0667;
    text(5, 0.7, [num2str(time) ' sec'], 'FontSize',10)
    drawnow();

    if recordvid == 1
        frame = getframe(vidfig);
        writeVideo(v,frame)
    end
end
% %  %%
% flpstart = 107;
% flpend = 116;
% d.fixedAxialBrightField(startidx+flpstart:startidx+flpend,:) = fliplr(d.fixedAxialBrightField(startidx+flpstart:startidx+flpend,:));
%%
gdata = d.autoAxialSignal(startidx:endidx, :);
bfdata = d.fixedAxialBrightField(startidx:endidx,:);
pad = NaN(size(gdata,1),5);

cat = [bfdata pad gdata];
figure()
imagesc(cat, [8000 20000])
colormap turbo
%%
ax = gca;
yt = round(ax.YTick/15);
for i = 1:length(yt)
    ylab(i) =  {num2str(yt(i))};
end

ax.YTickLabel = ylab;
ax.XTick = [130 390];
ax.XTickLabel = {'mCherry', 'GCaMP'};
ylabel('Time (s)');

exportgraphics(ax, strrep(savename, '.mp4', '.png'), 'Resolution', 300)


close(v)
%%
% wormdata = d;
% save(mfile, 'wormdata')
end