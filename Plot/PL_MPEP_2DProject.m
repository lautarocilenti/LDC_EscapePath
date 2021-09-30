function [iS] = PL_MPEP_2DProject(data)
%PL_MPEP 
theta = data.theta(data.minPhiIndex);
minPhi = data.phiSet{data.minPhiIndex}; 
minS = data.minS;
M = data.M;



[minS,iS] = min(minS);

t = minPhi{1};
phiAll = minPhi{2};
if M.plotFall
    tFall = minPhi{4};
    PhiFall = minPhi{5};
end
PL_Attractors(data.attractors);
hold on
for i = 1:M.dim/2
    phi = phiAll(:,(i-1)*2+1:(i)*2);

        color = [rand() rand() rand()];
        [tq,xq] = InterpolateToPhaseAngle(t,phi,M);
        plot(xq(:,1),xq(:,2),'o-','Color',color);
        hold on
        if M.plotFall
            phiFall = PhiFall(:,(i-1)*2+1:(i)*2);
            [tqFall,xqFall] = InterpolateToPhaseAngle(tFall,phiFall,M);
            plot([xq(end,1);xqFall(:,1)],[xq(end,2);xqFall(:,2)],'x--','Color',color);      
        end
    
end

axis([-4 4 -4 4])


% sprintf("MPEP with $\\nu$: %.2f, $\\theta$: %.3f, S: %.3f R: %.3e",M.Mrhs.nu,theta,minS,M.rIC)
%title(sprintf("MPEP with $\\nu$: %.2f, $\\theta$: %.3f, S: %.3f R: %.3e",M.Mrhs.nu,theta,minS,M.rIC),'interpreter','latex')
xlabel("$q_1$",'interpreter','latex')
ylabel("$q_2$",'interpreter','latex')

end

