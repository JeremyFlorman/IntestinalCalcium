function reCombine_Wormdata(filenames, outputdir, genotype)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pt = '';
for j = 1:length(filenames)
    pt = [pt '.'];
    disp(pt);

    wd = copyLoadClear(filenames{j}, 'C:\tmp');
    wormdata = wd.wormdata;
    clear('wd');
    

    fields = fieldnames(wormdata);
    for m = 1:length(fields)
        mergedstructure(j).(fields{m}) = wormdata.(fields{m}); % merge replicate data
    end


    if ~isfield(wormdata, 'autoAxialSignal') && isfield(wormdata, 'rawAxialSignal')
        mergedstructure(j).autoAxialSignal = autoFixSignal(wormdata.rawAxialSignal);
    end

    mergedstructure(j).filename = filenames(j);
    mergedstructure(j).genotype = genotype;
end


if isfield(mergedstructure,'bulkAboveBkg')
    mergedstructure = rmfield(mergedstructure, 'bulkAboveBkg');
end

if isfield(mergedstructure,'sumSignal')
    mergedstructure = rmfield(mergedstructure, 'sumSignal');
end

wormdata = mergedstructure;
structureSaveName = [outputdir '\' genotype '_mergedData.mat']
save(structureSaveName, 'wormdata');


disp("Done Combining wormdata")
