function [] = PL_L2MPEPPath(data,type,quantity)
%PL_L2PATH 


if nargin == 1
   type = "action"; 
   quantity = 1;
elseif nargin == 2
    quantity = 1;
end

%PL_MPEP 
if strcmp(type,"action")
    costSet = data.S;
else
    costSet = data.C;
end

[~,iC] = sort(costSet,'ascend');

M = data.M;
T = 2*pi/M.Mrhs.w;
FP = M.attractors;
vecnorm(FP,2,2)


for i = 1:quantity
    phi = data.phiSet{iC(i)};
    PL_PhiL2(phi,FP,T,M,true)
end


end

