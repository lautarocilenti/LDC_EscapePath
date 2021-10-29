function [iS] = PL_MPEP_2DProject(data)
%PL_MPEP 
theta = data.theta(data.minPhiIndex);
minPhi = data.phiSet{data.minPhiIndex}; 
minS = data.minS;
M = data.M;



[minS,iS] = min(minS);

 PL_Path_2DProject(data,iS);

end

