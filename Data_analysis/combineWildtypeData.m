function [] = combineWormdata(datafolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin <1
    %     datafolder = uigetdir('Z:\Calcium Imaging\Intestinal_Calcium_FreelyMoving')
    datafolder = 'Y:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\';
% datafolder = 'Y:\Calcium Imaging\Intestinal_Calcium\Rebekka\dec_mutants\';
end

 combinedDataFolder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\wildtype';


 
d = dir([datafolder '\**\*wildtype_*_wormdata.mat']);
nms = {d.name};
names = cell(1,length(nms));
dates = nan(1,length(nms));
for i = 1:length(nms)
    r = regexp(nms{i}, '_zfis178_', 'split');

   dates(i)=str2double(r{1});
    
    r = regexp(r{2}, '_', 'split');
    names(i) = {r(1)};
end
isremote = 1; 
genotypes = unique([names{:}]);
% genotypes = genotypes(2:end);

dateflag = dates>230410;

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
        

        if ~isfield(wormdata, 'include')
            mergedstructure(j).include = 1;
        end

        mergedstructure(j).filename = files2merge(j);
        mergedstructure(j).genotype = genotypes{k};
    end
    
    wormdata = mergedstructure;
    structureSaveName = [combinedDataFolder '\' genotypes{k} '_mergedData.mat']
    save(structureSaveName, 'wormdata');
    
    
    
    
    
end
