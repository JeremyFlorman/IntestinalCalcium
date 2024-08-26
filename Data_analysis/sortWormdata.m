% import wormdata file and create a variable called 'score' with the values
% to sort by.

% score = [];

l = score == 1;
m = score == 2;
h = score == 3; 

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
    for i = 1:length(data)
        datamean(i) = mean(data(i).autoAxialSignal,'all','omitnan');
    end

    [~, sortorder] = sort(datamean,'descend');
    for i = 1:length(sortorder)
    sortedData(i) = data(sortorder(i));
    end
end