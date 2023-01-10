function [autoAxialSignal] = autoFixSignal(axialSignal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

autoAxialSignal = axialSignal;
segment = ceil(size(autoAxialSignal, 2)*0.1);
for ii = 1:length(axialSignal)
    left = mean(axialSignal(ii,1:segment),'omitnan');
    right = mean(axialSignal(ii,end-segment:end),'omitnan');
    if left > right
        autoAxialSignal(ii,:) = fliplr(axialSignal(ii,:));
    end
end
end

