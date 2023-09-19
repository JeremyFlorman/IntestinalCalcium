folder  = 'X:\Calcium Imaging\Intestinal_Calcium\DMP_Mutants\';

d = dir(folder);
d = d(3:end);

for i = 22:length(d)
    datafolder = fullfile(d(i).folder,d(i).name);
    combineWormdata(datafolder)
end
