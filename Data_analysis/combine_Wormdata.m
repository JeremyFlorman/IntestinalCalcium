function combine_Wormdata(datadir,outputdir, settings)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

isremote = 1;
if nargin <1 % modify these paramaters for use without the app
    datadir = 'Y:\OAS\5-HT\+Food\';
    outputdir = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\Exogenous_Tyramine\receptor_mutants\wildtype-30mM-TA';

end

if isempty(outputdir)
    outputdir = datadir;
    isremote = 0;
end

% use this to combine genotypes that may be named slightly differently but
% should be combined. Enclose keywords in asterisks 
combineByKeyword =0;
keywords = {'*DA+Lactate*', '*DA-Lactate*', 'PhaC+Lactate', 'PhaC-Lactate'};



if combineByKeyword == 0

    dd = dir([datadir '\**\*_wormdata.mat' ]);

    names = cell(length(dd),1);
    dates = nan(length(dd),1);

    for i = 1:length(dd)
        tempname = split(dd(i).name, '_');
        names(i) = {tempname(end-2)};
    end

    genotypes = unique([names{:}]);

elseif combineByKeyword == 1
    genotypes = keywords;
end

%% add behavioral annotations from spreadsheet to wormdata
if settings.addBehaviorAnnotations == 1
    addBehaviorAnnotations(datadir, settings.framerate)
end



for i = 1:length(genotypes)
    d = dir([datadir '\**\*_' genotypes{i} '_*_wormdata.mat']);
    [datadir '\**\*_' genotypes{i} '_*_wormdata.mat']

    %% Exclude/include experiments based on date specified in dateflag
    nms = {d.name};
    dates = nan(length(nms),1);
    for j = 1:length(nms)
        r = regexp(nms{j}, '_zfis178_', 'split');
        dates(j)=str2double(r{1});
    end

    combineByDate = 0;
    if combineByDate == 1
    dateFlag = find(dates < 250101); %set date range to include in merged dataset. format is yymmdd
    d = d(dateFlag);
    end
    %% 

    disp(['Found ' num2str(length(d)) ' files with genotype ' genotypes{i}])

    mergedstructure = struct();

    pt = '';
    for j = 1:length(d)
        fileName = fullfile(d(j).folder, d(j).name);
        pt = [pt '.'];
        disp(pt);


        if isremote == 1
            wd = copyLoadClear(fileName, outputdir);
            wormdata = wd.wormdata;
            clear('wd');
        else
            load(fileName, 'wormdata');
        end


        fields = fieldnames(wormdata);
        for m = 1:length(fields)
            mergedstructure(j).(fields{m}) = wormdata.(fields{m}); % merge replicate data
        end


        if ~isfield(wormdata, 'autoAxialSignal') && isfield(wormdata, 'rawAxialSignal')
            mergedstructure(j).autoAxialSignal = autoFixSignal(wormdata.rawAxialSignal);
        end

        mergedstructure(j).filename = fileName;
        mergedstructure(j).genotype = genotypes{i};
    end


    if isfield(mergedstructure,'bulkAboveBkg')
        mergedstructure = rmfield(mergedstructure, 'bulkAboveBkg');
    end

    if isfield(mergedstructure,'sumSignal')
        mergedstructure = rmfield(mergedstructure, 'sumSignal');
    end

    wormdata = mergedstructure;
    structureSaveName = [outputdir '\' strrep(genotypes{i}, '*','') '_mergedData.mat']
    save(structureSaveName, 'wormdata');





end
disp("Done Combining wormdata")

% mergeControl(outputdir,controlname)
% disp("Done merging control")
%
% plot_MatchedControl(outputdir)
% disp('done plotting!!!!!')
