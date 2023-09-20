function [] = mergeControl(datafolder,controlname)
%Adds control data field to wormdata structure
%   this will take a control dataset specified by the text in 'controlname'
%   and adds it as a field to another dataset in the same folder. this is
%   assumed to be the mutant dataset.

if nargin <1
    datafolder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\5-HT\wildtype-30mM-5HT';
    controlname = 'wildtype-control';
end



flds = dir(datafolder);
dirs = [flds.isdir];


if nnz(dirs) == 2 % if we arent doing a loop just assign 1 to flds
    flds = 1;
else                    % if looping, remove hidden dirs.
    dirs(1:2) = 0;
    flds = flds(dirs);
    
    ctrlflag = nan(length(flds),1); 
    for q = 1:length(flds)
        ctrlflag(q) = ~strcmp(flds(q).name, controlname);
        
    end
    flds = flds(logical(ctrlflag));
 end

for i = 1:length(flds)

    if length(flds) == 1
        subfolder = datafolder
    else
        subfolder = fullfile(flds(i).folder, flds(i).name)
    end
    d = dir([subfolder '\**\*mergeddata.mat']);
    names = {d(:).name};

    cont = contains(names, controlname,"IgnoreCase",true);
    if nnz(cont) == 0
        disp('Cant find control genotype... check value of "controlname"')
    end

    mt = ~cont;
 
    mtstruct = load(fullfile(d(mt).folder, d(mt).name));
    wtstruct = load(fullfile(d(cont).folder, d(cont).name));

    controldata = wtstruct.wormdata;
    mtstruct.wormdata(1).controlData = controldata;

    wormdata = mtstruct.wormdata;
    savename = fullfile(d(mt).folder, d(mt).name)

    save(savename, "wormdata");
end

disp('Complete!')
end
