function [plotorder] = genotypeOrder(genotypes,parentfolder)
%genotypeOrder returns a vector containing the order of genotypes in a
%folder as specified in the genotypes input variable. 
%   Detailed explanation goes here
genotypes = {'wildtype', 'cca-1', 'egl-19(gf)','egl-19(lf)', 'trpa-1', ...
    'flr-1', 'eat-2', 'tph-1', 'inx-2', 'inx-16', 'itr-1','unc-43(n498sd)',...
    'unc-43(sa200)', 'unc-43(e408)', 'itr-1-unc-43','dec-1','dec-7','dec-9','dec-10'};

parentfolder = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants';
dd = dir(parentfolder);
ddirflag = [dd.isdir];
dd = dd(ddirflag);
dd = dd(3:end);

plotorder = nan(length(genotypes),1);

for i = 1:length(genotypes)
  plotorder(i) = find(strcmp({dd.name},genotypes{i}));


end

% for i = 1:length(genotypes)
% contains([dd.name])
% end

end