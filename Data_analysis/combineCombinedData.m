datafolder = 'C:\Users\Jeremy\Dropbox\Intestinal Calcium Paper\Data\Figure 2';
combinedDataFolder = 'C:\Users\Jeremy\Dropbox\Intestinal Calcium Paper\Data\Figure 2\wildtype';

isremote = 0;


d = dir([datafolder '\**\*wildtype_mergedData.mat']);
nms = {d.name};

genotypes = {'wildtype'};
ii = 1;

for k = 1:length(genotypes)
    mergedstructure = struct();
    replicates = find(contains(nms, genotypes{k})); % index of replicates for a given genotype
    files2merge = cell(length(replicates),1);       % cell array to hold path for replicates
    disp(genotypes{k})
    disp(nms(contains(nms, genotypes{k})))



    for j = 1:length(replicates)
        files2merge(j) = {fullfile(d(replicates(j)).folder, d(replicates(j)).name)};

        tic

        if isremote == 1
            wd = copyLoadClear(files2merge{j}, combinedDataFolder);
            wormdata = wd.wormdata;
            clear('wd');
        else
            load(files2merge{j});
        end

        toc

        fields = fieldnames(wormdata);

        for expi = 1:length(wormdata)

            for m = 1:length(fields)
                mergedstructure(ii).(fields{m}) = wormdata(expi).(fields{m}); % merge replicate data
            end


            if ~isfield(wormdata, 'autoAxialSignal') && isfield(wormdata, 'rawAxialSignal')
                mergedstructure(ii).autoAxialSignal = autoFixSignal(wormdata.rawAxialSignal);
            end


            if ~isfield(wormdata, 'include')
                mergedstructure(ii).include = 1;
            end
            ii = ii+1;
        end

    end

    wormdata = mergedstructure;
    structureSaveName = [combinedDataFolder '\' genotypes{k} '_mergedData.mat']
    save(structureSaveName, 'wormdata');





end
