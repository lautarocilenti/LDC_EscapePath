function [] = PL_L2MPEPPathMultipleCandidates(data)
%PL_L2PATH Summary of this function goes here
%   Detailed explanation goes here
M = data.M;


[M] = GetFixedPoints(M);

s = data.S;
theta = data.theta;



nMinPhi = 3;
L = 1;


for i = 1:nMinPhi
    [~,iSort] = sort(s,'ascend');
    minPhiIndex(i) = iSort(1);
    distance = vecnorm(theta-theta(:,minPhiIndex(i)),2,1);
    iInRange = distance <=L;
    s(iInRange) = s(iInRange)+1000;
end




FP = M.Mrhs.FixedPoints.FP;
t = M.tspan;
T = 2*pi/M.Mrhs.w;

data.S(minPhiIndex)

for i = 1:nMinPhi
    phi = data.phiSet{minPhiIndex(i)};
    PL_PhiL2(phi,FP,T,M,data.fullpath)
    hold on

end
title(sprintf("%.8f,%.8f,%.8f",data.S(minPhiIndex)))


end
