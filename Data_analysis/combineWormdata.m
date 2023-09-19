function [] = combineWormdata(datafolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <1
        datafolder = 'Y:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\Receptor_Mutants\wildtype-30mM-TA'
end
 
folderparts = regexp(datafolder, '\','split');

combinedDataDir = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\Exogenous_Tyramine\receptor_mutants\';

combinedDataFolder = [combinedDataDir, folderparts{end}]
controlname = 'wildtype-control';

if ~isfolder(combinedDataFolder)
    mkdir(combinedDataFolder) 
end

 
d = dir([datafolder '\**\*wormdata.mat']);

%% Exclude/include experiments based on date specified in dateflag
nms = {d.name};
dates = nan(1,length(nms));
for i = 1:length(nms)
    r = regexp(nms{i}, '_zfis178_', 'split');

    dates(i)=str2double(r{1});
end

dateflag = dates>210410; %set date range to include in merged dataset. format is yymmdd

d = d(dateflag);

%%

nms = {d.name};
names = cell(1,length(nms));
for i = 1:length(nms)
    r = regexp(nms{i}, 'zfis178_', 'split');
    
    r = regexp(r{2}, '_', 'split');
    names(i) = {r(1)};
end
isremote = 1;
genotypes = unique([names{:}]);

if length(genotypes)>1 && strcmp(genotypes{2}, 'itr-1-unc-43')
     genotypes = genotypes(2:end);
end

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
        for m = 1:length(fields)
            mergedstructure(j).(fields{m}) = wormdata.(fields{m}); % merge replicate data
        end
        
        
        if ~isfield(wormdata, 'autoAxialSignal') && isfield(wormdata, 'rawAxialSignal')
            mergedstructure(j).autoAxialSignal = autoFixSignal(wormdata.rawAxialSignal);
        end
        
        mergedstructure(j).filename = files2merge(j);
        mergedstructure(j).genotype = genotypes{k};
    end
    
    wormdata = mergedstructure;
    structureSaveName = [combinedDataFolder '\' genotypes{k} '_mergedData.mat']
    save(structureSaveName, 'wormdata');
    
    
    
    
    
end
disp("Done Combining wormdata")

mergeControl(combinedDataFolder,controlname)
disp("Done merging control")

plot_MatchedControl(combinedDataFolder)
disp('done plotting!!!!!')
