function [iS] = PL_MPEP(data)
%PL_MPEP 
theta = data.theta(data.minPhiIndex);
minPhi = data.phiSet{data.minPhiIndex}; 
minS = data.minS;
M = data.M;
rng(2)
for k = 1:2
    
    %temp additions
    if k == 1
    ii = find(data.theta<2);
    else
    ii = find(data.theta>2);   
    end
    theta = data.theta(ii);
    phiSet = data.phiSet(ii); 
    [minS,iPhi] = min(data.S(ii));
    theta = theta(iPhi)
    minPhi = phiSet{iPhi};
    %end temp
    phi = minPhi;
    color = [rand() rand() rand()];

    [minS,iS] = min(minS);
    t = minPhi{1};
    phi = minPhi{2};

    if strcmp(M.rhsString,'Unforced');
        plot(phi(:,1),phi(:,2),'linewidth',2)
    else strcmp(M.rhsString,'Duffing');
        color = [rand() rand() rand()];
        [tq,xq] = InterpolateToPhaseAngle(t,phi,M);
        if k == 1
        plot(xq(:,1),xq(:,2),'o-','Color',color);
        else
            plot(xq(:,1),xq(:,2),'x-','Color',color);
        end
%         if M.plotFall
%             tFall = minPhi{4};
%             phiFall = minPhi{5};
%             [tqFall,xqFall] = InterpolateToPhaseAngle(tFall,phiFall,M);
%             plot([xq(end,1);xqFall(:,1)],[xq(end,2);xqFall(:,2)],'x--','Color',color);      
%         end
    end
    hold on
end
title(sprintf("MPEP with $\\nu$: %.2f, $\\theta$: %.3f, S: %.3f R: %.3e",M.Mrhs.nu,theta,minS,M.rIC),'interpreter','latex')
xlabel("$q_1$",'interpreter','latex')
ylabel("$q_2$",'interpreter','latex')

end

