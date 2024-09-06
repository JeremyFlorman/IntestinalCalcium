% import wormdata file and create a variable called 'score' with the values
% to sort by.

% score = [];
d = dir('C:\Users\Jeremy\Dropbox\ins-3 intestinal calcium fig\data\pooledData\*mergedData.mat');
excelpath = "Z:\Calcium Imaging\Intestinal_Calcium\Exogenous_Tyramine\TA Ca2+ Response.xlsx";

scores = readmatrix(excelpath);
s(1) = {scores(:,28)}; %itr-1
s(2) ={scores(:,7)}; %lgc-55
s(3) ={scores(:,38)}; %quad
s(4) ={scores(:,24)}; %tyra-3
s(5) ={scores(:,3)}; %wt

for k = 1:length(d)
    datapath = fullfile(d(k).folder,d(k).name);
    load(datapath)
    
    score = s{k};
    nanflag = isnan(score);
    score = score(~nanflag);

    sortedData = sortData(wormdata, score);
    wormdata = sortedData;
    save(datapath,"wormdata")
end


function sortedData=sortData(wormdata,score)
l = score == 1;
m = score == 2; 
h = score ==3;
 
low = wormdata(l);
med = wormdata(m);
high = wormdata(h);

%%
clear('sortedData')

sortedLow = sortByMean(low);
sortedMed = sortByMean(med);
sortedHigh = sortByMean(high);

for i = 1:length(sortedHigh)
    sortedData(i) = sortedHigh(i);
end

buffer1 = length(sortedHigh);
for i = 1:length(sortedMed)
    sortedData(i+buffer1) = sortedMed(i);
end

buffer2 = length(sortedHigh)+length(sortedMed);
for i = 1:length(sortedLow)
    sortedData(i+buffer2) = sortedLow(i);
end


    function sortedData = sortByMean(data)
        datamean = nan(length(data),1);
        for k = 1:length(data)
            datamean(k) = mean(data(k).autoAxialSignal,'all','omitnan');
        end

        [~, sortorder] = sort(datamean,'descend');
        for k = 1:length(sortorder)
            sortedData(k) = data(sortorder(k));
        end
    end
end