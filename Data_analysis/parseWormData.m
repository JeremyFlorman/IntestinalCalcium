function [mtdata, wtdata, settings] = parseWormData(datapath,settings)
%parseWormData - Returns mutant and control data from a mergedData
%file.
%   if datapath is a file location, that file will be loaded. if datapath
%   is a structure that has already been loaded, that structure will be
%   used instead of loading. The "mutant" is stored in the structure
%   wormData. Matching control data is extracted from
%   wormData(1).controlData. Use mergeControl.m to store control data.

% datapath = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\DMP_mutants\egl-19(lf)\egl-19(lf)_mergedData.mat';

if nargin<2
settings = returnPlotSettings();
end

normalize = settings.normalize;

if ischar(datapath) || isstring(datapath) 
    mtdata = load(datapath);
    mtdata = mtdata.wormdata;
elseif isstruct(datapath)
    mtdata = datapath;
end

if isfield(mtdata, 'include')
    mtdata = mtdata(logical([mtdata.include]));
end


if isfield(mtdata, 'controlData')
    wtdata = mtdata(1).controlData;
    if isfield(wtdata, 'include')
        wtdata = wtdata(logical([wtdata.include]));
    end



    % normalize bulk signal by dividing by the mean wildtype bulk signal.
    if normalize == 1
        wtbulksig = cell2mat({wtdata(:).bulkSignal});
        wtmedian = median(wtbulksig,'all','omitnan');
        for i = 1:length(mtdata)
            mtdata(i).normalizedSignal = mtdata(i).bulkSignal./wtmedian-1;
        end

        for i = 1:length(wtdata)
            wtdata(i).normalizedSignal = wtdata(i).bulkSignal./wtmedian-1;
        end
    end

else
    if normalize == 1
        mtbulksig = cell2mat({mtdata(:).bulkSignal});
        mtmedian = median(mtbulksig,'all','omitnan');
        for i = 1:length(mtdata)
            mtdata(i).normalizedSignal = mtdata(i).bulkSignal./mtmedian-1;
        end
    end
    wtdata = mtdata;
end





%%
    function [normalizedSignal] = normalizeSignal(signal2Normalize, controlSignal)
        avgControlSignal = median(cell2mat({controlSignal(:).bulkSignal}),'all','omitnan');
        for  i2 = 1:length(signal2Normalize)
        signal2Normalize(i2).normalizedSignal = signal2Normalize(i2).bulkSignal./avgControlSignal;
        end

        
    end

        
end