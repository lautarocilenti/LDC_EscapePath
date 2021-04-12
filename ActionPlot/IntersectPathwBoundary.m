function [iIntersect] = IntersectPathwBoundary(basinInterpolant,phiSet,M)
%INTERSECTPATHWBOUNDARY
iIntersect = -1*ones(M.nIC,1);
for iPhi = 1:M.nIC
    q = phiSet(:,1:length(M.qo),iPhi);
    try 
        iIntersect(iPhi) = min(find(round(basinInterpolant(q),0) == 2));
    catch
%         warning('Solution for iPhi did not intersect boundary');
         error('Solution for iPhi did not intersect boundary');
    end
end


end

