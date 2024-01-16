function combine_Wormdata(datadir,outputdir, controlname)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

isremote = 1;
if nargin <1 % modify these paramaters for 
    datadir = 'Y:\OAS\5-HT\+Food\';
    outputdir = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\Exogenous_Tyramine\receptor_mutants\wildtype-30mM-TA';
    
end



if isempty(outputdir)
    outputdir = datadir;
    isremote = 0;
end



% genotypes = {mtname controlname};

dd = dir([datadir '\**\*_wormdata.mat' ]);

names = cell(length(dd),1);

for i = 1:length(dd)
    tempname = split(dd(i).name, '_');
    names(i) = {tempname(end-2)};
end

genotypes = unique([names{:}]);






for i = 1:length(genotypes)  
    d = dir([datadir '\**\*' genotypes{i} '_*_wormdata.mat']);

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
    
    wormdata = mergedstructure;
    structureSaveName = [outputdir '\' genotypes{i} '_mergedData.mat']
    save(structureSaveName, 'wormdata');
    
    
    
    
    
end
disp("Done Combining wormdata")

% mergeControl(outputdir,controlname)
% disp("Done merging control")
% 
% plot_MatchedControl(outputdir)
% disp('done plotting!!!!!')
