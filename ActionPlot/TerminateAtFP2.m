function [value, isterminal, direction] = TerminateAtFP(t, y, Mrhs)
value = [1;1;1];

if mod(Mrhs.w*t,2*pi) <= Mrhs.psiEps %once a period
    Mrhs.temp.told = t;
    iFP = find(xInSpheres(Mrhs.FixedPoints.FP,Mrhs.rA,y'));
    LiFP = length(iFP);
    if LiFP == 1 %inside of sphere of one fixed point
        if Mrhs.FixedPoints.Stability(iFP) == 1 %attractor
            if xInSpheres(Mrhs.FixedPoints.FP(Mrhs.iA,:),Mrhs.rA,y') %initial attractor
%                 status = "InitialBasin"
                value(1) = 0;
            else
%                 status = "AlternateBasin"
                value(2) = 0;
            end
        else %saddle
            if xInSpheres(Mrhs.FixedPoints.FP(iFP,:),Mrhs.rS,y')
%                 status = "Saddle";
                value(3) = 0;
            else
%                 status = "Nothing";
            end
        end
    elseif length(iFP) > 1
        error("Multiple FP near solution")
    else
%         status = "Nothing";
    end
else
    %nothing
end

isterminal = [1;1;1];   % Stop the integration
direction  = [0;0;0];
end

function [res] = xInSpheres(C,R,x)
    res = vecnorm(C-x,2,2)<=R;
end