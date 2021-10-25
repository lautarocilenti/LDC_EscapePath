function [iBestNewValue] = BestNewSteps(sNewSearch,iNewSearch)
%EVALUATENEWSTEP 
    iNewSearchUnique = unique(iNewSearch);
    iBestNewValue = zeros(1,length(iNewSearchUnique));
    for i = 1:length(iNewSearchUnique)
        ii = find(iNewSearchUnique(i) == iNewSearch);
        [~,imin] = min(sNewSearch(ii));
        iBestNewValue(i) = ii(imin);
    end

end

