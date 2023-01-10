[name, folder]  = uigetfile('C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants');

toggleMutant = [5 6 7 8 9];
toggleControl = [];

 wormdata = load(fullfile(folder, name)); 
 wormdata = wormdata.wormdata; 
for i = 1:length(toggleMutant)
 wormdata(toggleMutant(i)).include = abs(wormdata(toggleMutant(i)).include-1);
end

 if ~isempty(toggleControl) 
     wormdata(1).controlData(toggleControl).include = abs(wormdata(1).controlData(toggleControl).include-1);
 end

save(fullfile(folder, name), 'wormdata')

% 
% mf = dir([folder '\**\*mergedData.mat']);
% 
% for i = 1:length(mf)
%    
%     wormdata =wormdata.wormdata;
%         for k =1:size(wormdata,2)
%             wormdata(k).include = 1;
%         end
%         save(fullfile(mf(i).folder, mf(i).name), 'wormdata')
% end