folder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\Exogenous_Tyramine\receptor_mutants\';

mf = dir([folder '\**\*mergedData.mat']);

for i = 1:length(mf)
    load(fullfile(mf(i).folder, mf(i).name));
    
        for k =1:size(wormdata,2)
            wormdata(k).include = 1;
        end
        save(fullfile(mf(i).folder, mf(i).name), 'wormdata')
end




