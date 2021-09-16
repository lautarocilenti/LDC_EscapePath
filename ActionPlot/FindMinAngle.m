function [minAngle] = FindMinAngle(S,phiSet,M)
%FINDMINANGLE Summary of this function goes here
%   Detailed explanation goes here
ii = find(S>0);
[~,iPhi] = min(S(ii));
minAngle = M.theta(iPhi);
end

