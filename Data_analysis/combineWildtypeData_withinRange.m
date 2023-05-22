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
dates = nan(1,length(nms));
for i = 1:length(nms)
    r = regexp(nms{i}, '_zfis178_', 'split');

    dates(i)=str2double(r{1});
end

dateflag = dates>230410; %set date range to include in merged dataset. format is yymmdd

d = d(dateflag);

isremote = 1;
genotype = 'wildtype';
mergedstructure = struct();
files2merge = cell(length(d),1);       % cell array to hold path for replicates



for j = 1:length(d)
    files2merge(j) = {fullfile(d(j).folder, d(j).name)};

    if isremote == 1
        wd = copyLoadClear(files2merge{j}, combinedDataFolder);
        wormdata = wd.wormdata;
        clear('wd');
    else
        load(files2merge{j});
    end



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
    mergedstructure(j).genotype = genotype;
end

wormdata = mergedstructure;
structureSaveName = [combinedDataFolder '\' genotype '_mergedData.mat']
save(structureSaveName, 'wormdata');






