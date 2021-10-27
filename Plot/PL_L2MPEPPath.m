function [] = PL_L2MPEPPath(data)
%PL_L2PATH Summary of this function goes here
%   Detailed explanation goes here
M = data.M;


[M] = GetFixedPoints(M);
theta = data.theta(data.minPhiIndex);
minPhi = data.phiSet{data.minPhiIndex}; 
minS = data.minS;

FP = M.Mrhs.FixedPoints.FP;
t = M.tspan;
T = 2*pi/M.Mrhs.w;

phi = minPhi;

PL_PhiL2(phi,FP,T,M,data.fullpath)



end

