function [] = PL_L2MPEPPath(data)
%PL_L2PATH Summary of this function goes here
%   Detailed explanation goes here
M = data.M;


  [M] = GetFixedPoints(M); 
s = data.S;
[sSorted,iSortS] = sort(s,'ascend');

phiSet = data.phiSet(iSortS);

FP = M.Mrhs.FixedPoints.FP;
t = M.tspan;
T = 2*pi/M.Mrhs.w;
for i = 1:100:length(phiSet)
    phi = phiSet{i};

    PL_PhiL2(phi,FP,T,M)
    hold off
end






end

