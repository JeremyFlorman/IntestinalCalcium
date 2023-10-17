% [file, path] = uigetfile('C:\Users\Jeremy\Dropbox\Intestinal Calcium Paper\Data\*mergedData.mat');
% load(fullfile(path,file));
load("C:\Users\Jeremy\Dropbox\Intestinal Calcium Paper\Data\Figure 2\wildtype\wildtype_mergedData.mat");


figure();
t = tiledlayout(2,2);
bulkAx = nexttile([1 2]);
axAx = nexttile([1 2]);

includeVector = nan(length(wormdata),1);


for i = 1:length(wormdata)
    bulkSignal = fillmissing(wormdata(i).bulkSignal-wormdata(i).backgroundSignal, 'movmedian',100);

    plot(bulkSignal, 'Parent', bulkAx);
    bulkAx.XLim = [0 length(bulkSignal)];
    bulkAx.YLim = [-500 8000];


    backgroundMatrix = repmat(wormdata(i).backgroundSignal,1,size(wormdata(i).autoAxialSignal,2));
    axsig = wormdata(i).autoAxialSignal-backgroundMatrix;

    axylimit = [0 35000];
    imagesc(smoothdata(axsig,'movmedian',60)', axylimit)
    axAx.XLim = [0 length(axsig)];
    axAx.YDir = "reverse";

    [~, name, ~] = fileparts(wormdata(i).filename);
    title(t, strrep(name, '_', ' '));

    txt = input('accept experiment? Y/N, type exit to quit', 's' );

    if strcmpi(txt, 'y')
        includeVector(i) = 1;
    elseif strcmpi(txt, 'n')
        includeVector(i) = 0;
    elseif strcmpi(txt, 'exit')
        break
    end

end
%%
for j = 1:length(includeVector)
    wormdata(j).include = includeVector(j);
end
