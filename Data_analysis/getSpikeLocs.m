function [spikeProperties] = getSpikeLocs(datapath, getControl,settings)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


% peakthreshold = 500;
% datapath = 'C:\Users\Jeremy\Desktop\Calcium Imaging\FreelyMoving_Data\combinedData\Exogenous_Tyramine\receptor_mutants\lgc-55-30mM-TA\lgc-55-30mM-TA_mergedData.mat';
% getControl =0;


[mtdata, wtdata, settings] = parseWormData(datapath);

peakdistance = settings.peakdistance;
peakthreshold = settings.peakthreshold;
peakwidth = settings.peakwidth;
normalize = settings.normalize;
sortType = settings.sortType;
sortDir = settings.sortDir;

mtamp = [];
mtlocs = [];
wtamp = [];
wtlocs = [];
mtsignal = {};
wtsignal = {};
mtinterval = [];
wtinterval = [];

mtnum = nan(length(mtdata),1);
wtnum = nan(length(wtdata),1);
mtAvAmp = nan(length(mtdata),1);
wtAvAmp = nan(length(wtdata),1);

mtsort = nan(length(mtdata),1);
wtsort = nan(length(wtdata),1);

% if normalize == 1
% idx = ~cellfun(@isempty,{mtdata.normalizedSignal});
% mtbulksig = mtdata(idx).normalizedSignal;
% end

for i = 1:length(mtdata)

    if normalize == 1
        sig = mtdata(i).normalizedSignal;
    else
        sig = mtdata(i).bulkSignal;
    end


    terpdata = fillmissing(sig, 'movmedian',100);
    [tempamp, templocs] = findpeaks(terpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance,'MinPeakWidth',peakwidth);

    tempint = [];
    if length(templocs)>1
        tempint = diff(templocs)/15;
    end

    mtinterval = vertcat(mtinterval, tempint);
    mtamp = vertcat(mtamp, tempamp);
    mtlocs = vertcat(mtlocs, templocs);
    mtsignal(i) = {terpdata};
    if ~isempty(templocs)
        mtnum(i) = length(templocs);
        mtAvAmp(i) = mean(tempamp);
    else
        mtnum(i) = 0;
        mtAvAmp(i) = 0;
    end


end

if getControl == 1

    for i = 1:length(wtdata)

        if normalize == 1
            sig = wtdata(i).normalizedSignal;
        else
            sig = wtdata(i).bulkSignal;
        end

        terpdata = fillmissing(sig, 'movmedian',150);
        [tempamp, templocs] = findpeaks(terpdata, 'MinPeakProminence', peakthreshold, 'MinPeakDistance',peakdistance,'MinPeakWidth',peakwidth);

        tempint = [];
        if length(templocs)>1
            tempint = diff(templocs)/15;
        end
        

        if i == 1 
        wtinterval = tempint;
        wtamp = tempamp;
        wtlocs = templocs;
        elseif i>1
        wtinterval = vertcat(wtinterval, tempint);
        wtamp = vertcat(wtamp, tempamp);
        wtlocs = vertcat(wtlocs, templocs);
        end

        wtsignal(i) = {terpdata};

        if ~isempty(templocs)
            wtnum(i) = length(templocs);
            wtAvAmp(i) = mean(tempamp,'omitnan');
        else
            wtnum(i) = 0;
            wtAvAmp(i) = 0;
        end
    end
end



if sortType == 0                % dont sort
    mtsort = 1:length(mtdata);
elseif sortType == 1            % sort by number of spikes
    mtsorttype = mtnum;
    [~, mtsort] = sort(mtsorttype,sortDir);
elseif sortType == 2            % sort by spike amplitude
    mtsorttype = mtAvAmp;
    [~, mtsort] = sort(mtsorttype,sortDir);
end


spikeProperties.mtsort = mtsort;
spikeProperties.mtAmp = mtamp;
spikeProperties.mtLocs = mtlocs;
spikeProperties.mtSignal = [mtsignal{:}];
spikeProperties.mtInterval = mtinterval;

if getControl == 1


    if sortType == 0                % dont sort
        wtsort = 1:length(wtdata);
    elseif sortType == 1            % sort by number of spikes
        wtsorttype = wtnum;
        [~, wtsort] = sort(wtsorttype,sortDir);
    elseif sortType == 2            % sort by spike amplitude
        wtsorttype = wtAvAmp;
        [~, wtsort] = sort(wtsorttype,sortDir);
    end

    
    spikeProperties.wtsort = wtsort;
    spikeProperties.wtAmp = wtamp;
    spikeProperties.wtLocs = wtlocs;
    spikeProperties.wtSignal = wtsignal;
    spikeProperties.wtInterval = wtinterval;
end

end