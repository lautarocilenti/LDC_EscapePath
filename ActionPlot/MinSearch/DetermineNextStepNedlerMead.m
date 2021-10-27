function [thetaNewNonGradient,iNonGradientSearch,stepDirection] = DetermineNextStepNedlerMead(thetaCurrent,sCurrent,iCurrent,stepSize,State,M)
%DERTERMINENEXTSTEPNELDERMEAD 
%simplex state
% 1 - perform standard reflection
% 2 - Best last step, continue along past direction
% 3 - outside contract the simplex
% 4 - inside contract the simplex
% 5 - shrink the simplex
% 6 - do nothing


n = M.dim+1;
d = M.dim;
L = size(thetaCurrent,2);
Ln = L/n;
if M.xcoordinates | floor(Ln) ~= Ln
    error("Dimension issue in Nelder Mead algorithm")
end

simplexState = State.simplexState;

c = 1;
for i = 1:Ln
    j1 = (i-1)*n+1; j2 = (i)*n;
    simplex = thetaCurrent(:,j1:j2);

    simplexS = sCurrent(j1:j2);
    simplexI = iCurrent(j1:j2);
    [~,is] = sort(simplexS,"descend");
    iSearch = is(1);
    simplexCentroid = mean(simplex(:,is(2:end)),2);
    iBest = is(end);
    dir = simplexCentroid-simplex(:,iSearch);
    

    if simplexState(i) == 1;
        dir = dir; %standard reflection
    elseif simplexState(i) == 2;
        dir = 2*dir; %double reflection
    elseif simplexState(i) == 3;
        dir = .5*dir; %outside contraction
    elseif simplexState(i) == 4;
        dir = -.5*dir; %inside contraction
    elseif simplexState(i) == 5;
        %shrink
        for k = 1:n
            if k~=iBest
                 thetaNewNonGradient(:,c) = .5*(simplex(:,k)+simplex(:,iBest));
                 iNonGradientSearch(c) = simplexI(k);
                 stepDirection(:,c) = .5*(simplex(:,iBest)- simplex(:,k));
                 c = c+1;
            end
        end
        continue
    else
        continue %no more steps
    end
        
    
    thetaNewNonGradient(:,c) = simplexCentroid+dir;
    iNonGradientSearch(c) = simplexI(iSearch);
    stepDirection(:,c) = dir;
    c = c+1;
    
end


end

