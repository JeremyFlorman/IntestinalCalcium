function [patchX, patchY] = shadedFoodPatches(wormdata, ylims)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin<1
    wormdata = evalin("caller",'wormdata');
    ylims = [0 0.4169];
end

if isfield(wormdata, 'onFood')

    onFood = wormdata.onFood;
    offFood = wormdata.offFood;

    if length(offFood)<length(onFood)
        offFood(end+1) = length(wormdata.bulkSignal);
    end

    patchX = nan(4,length(onFood));
    patchY = nan(4,length(onFood));


    for i =1:length(onFood)
        patchX(1:2, i) = repmat(onFood(i),2,1); % left hand x coords
        patchY(1:2, i) = [ylims(1); ylims(2)];  % left hand y coords

        patchX(3:4,i) = repmat(offFood(i),2,1); % right hand x coords
        patchY(3:4, i) = [ylims(2); ylims(1)];                  % left hand y coords
    end

    % p = patch(patchX/15/60, patchY, [0.93 0.69 0.13], 'FaceAlpha', 0.5);

end
end