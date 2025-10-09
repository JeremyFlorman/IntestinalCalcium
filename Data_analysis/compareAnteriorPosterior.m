data = wildtypenoFood3hrData;

dataSize = size(data(1).autoAxialSignal);
anteriorSignal = nan(dataSize(1),length(data));
posteriorSignal = nan(dataSize(1),length(data));

timePre = 60; % seconds before alignment event
antStart = 20;
antEnd = antStart+dataSize(2)*0.1;

postStart = dataSize(2)- dataSize(2)*0.1;
postEnd = dataSize(2);

for i = 1:length(data)

    antMean = mean(data(i).autoAxialSignal(:,antStart:antEnd),2, 'omitmissing');
    postMean = mean(data(i).autoAxialSignal(:,postStart:postEnd),2, 'omitmissing');
    % imagesc(data(i).autoAxialSignal(:,antStart:antEnd))
    [antR, ~] = find(~isnan(antMean), 1, "first");
    antFZero = mean(antMean(antR:antR+15),1);

    [postR, ~] = find(~isnan(postMean), 1, "first");
    postFZero = mean(postMean(postR:postR+15),1);

    antDeltaF = (antMean-antFZero)/antFZero;
    postDeltaF = (postMean-postFZero)/postFZero;

    antZScore = (antDeltaF-mean(antDeltaF, 'omitmissing'))./std(antDeltaF, 'omitmissing');
    postZScore = (postDeltaF-mean(postDeltaF, 'omitmissing'))./std(postDeltaF, 'omitmissing');
    
    anteriorSignal(1:dataSize(1),i) = antMean;
    posteriorSignal(1:dataSize(1),i) = postMean;
end



time = linspace(-timePre,(dataSize(1)/15)-timePre, dataSize(1))';