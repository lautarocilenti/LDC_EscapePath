function [S] = ExtractAction(phiSet)
%EXTRACTACTION 
S = zeros(1,length(phiSet));
for iPhi = 1:length(phiSet)
   S(1,iPhi) = phiSet{iPhi}{6}; 
end

end

