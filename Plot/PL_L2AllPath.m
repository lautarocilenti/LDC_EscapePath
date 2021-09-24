function [] = PL_L2AllPath(data)
%PL_L2PATH Summary of this function goes here
%   Detailed explanation goes here
M = data.M;


  [M] = GetFixedPoints(M);
% theta = data.theta;
PhiSet = data.phiSet; 
% minS = data.minS;

FP = M.Mrhs.FixedPoints.FP;
t = M.tspan;
T = 2*pi/M.Mrhs.w;


for i = 1:100:length(PhiSet)
    phi = PhiSet{i};
    PL_PhiL2(phi,FP,T,M);
end



end

