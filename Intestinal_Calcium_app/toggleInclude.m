file = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\egl-19(gf)\egl-19(gf)_mergedData.mat';
[folder, name, ~] = fileparts(file);

% [name, folder]  = uigetfile('C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\Sammy\tir-1(qd4)\');

toggleMutant = [1 2 5 6];
toggleControl = [1 2 5 6];
 
 wormdata = load(fullfile(folder, name)); 
 wormdata = wormdata.wormdata; 
for i = 1:length(toggleMutant)
 wormdata(toggleMutant(i)).include = abs(wormdata(toggleMutant(i)).include-1);
end

 if ~isempty(toggleControl) 
     firstIncluded = find([wormdata.include],1);
     controlIndex = [];
     for k = 1:length(wormdata)
         if ~isempty(wormdata(k).controlData)
             controlIndex = k;
         end
     end
                
     
     controlData = wormdata(controlIndex).controlData;
     wormdata(controlIndex).controlData = [];

     for i = 1:length(toggleControl)
         
         controlData(toggleControl(i)).include = abs(controlData(toggleControl(i)).include-1);
     end
     wormdata(firstIncluded).controlData = controlData;

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