function [S] = MaximumTime(phiSet,M)
%MaximumTime 

S = zeros(1,length(phiSet));

for iPhi = 1:length(phiSet)
   t = phiSet{iPhi}{1};
   S(1,iPhi) = -max(t);
end


end

