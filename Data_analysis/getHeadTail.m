function [overlapX,overlapY,lineXY] = getHeadTail(xSegPts, ySegPts,xEvalPts, outline)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

                [p, ~, mu] =polyfit(xSegPts,ySegPts,1);  %fit those segments
                f = polyval(p,xEvalPts, [], mu);            %Evaluate that fit
                    
                [yOutPoints, xOutPoints] = find(outline);  %Get rotated outline points to compare
                linex = ceil(xEvalPts);
                liney = ceil(f);

                outlineXY = [xOutPoints yOutPoints];
                lineXY = [linex; liney]';

                 overlap = ismember(outlineXY,lineXY,'rows');
                 overlapX = outlineXY(overlap,1);
                 overlapY = outlineXY(overlap,2);
end