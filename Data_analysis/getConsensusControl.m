mtdir = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\OAS\5-HT\workingFigures\Fig3 - pumping mutants';
controlname = 'wildtype+5HT';

d = dir([mtdir '\**\*mergedData.mat']);
n = {d.name};
wtflag = contains(n,controlname, 'IgnoreCase', true);

d= d(~wtflag);
mergedControl = struct();

for i = 1:length(d)
    load(fullfile(d(i).folder,d(i).name));
    if isfield(wormdata, 'controlData')
        idx = find(~cellfun(@isempty,{wormdata.controlData}));
        tempcontrol = wormdata(idx).controlData;
        fields = fieldnames(tempcontrol);
        if i == 1
            mergedControl = tempcontrol;
        else
            mIdx = length(mergedControl);
            for j = 1:length(tempcontrol)
                for k = 1:length(fields)
                    aField = fields{k};
                    mergedControl(mIdx+j).(aField) = tempcontrol(j).(aField);
                end
            end
        end
    end
end

wormdata = mergedControl;
savedir = [mtdir '\wildtype+5HT'];
mkdir(savedir)
save([savedir '\wildtype+5HT_mergedData.mat'], "wormdata");
