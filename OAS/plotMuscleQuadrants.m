foldername = 'Y:\OAS\myo-3GCaMP\231204_zfex813_wildtype-10x+Tap_4\231204_zfex813_wildtype-10x+Tap_4_wormdata.mat';

load(foldername);
%%
topsignal = wormdata.topAxialSignal;
botsignal = wormdata.bottomAxialSignal;
nSegments = 12;

topTraces = nan(length(topsignal),nSegments);
botTraces = nan(length(botsignal),nSegments);
segmentEdges = floor(linspace(1,size(topsignal,2), nSegments+1));
time = linspace(0,length(topsignal)/15/60,length(topsignal))';



for i = 1:nSegments
    padVal = 75;
    padding = nSegments*padVal-i*padVal;
    topTraces(:,i) = mean(topsignal(:,segmentEdges(i):segmentEdges(i+1)),2)+padding;
    botTraces(:,i) = mean(botsignal(:,segmentEdges(i):segmentEdges(i+1)),2)+padding;
end


imagesc(topsignal)
plot(time,fliplr(topTraces))

ax = gca;
ax.ColorOrder = tab10(nSegments);