function [iS] = PL_Path_2DProject(data,index)
%PL_MPEP 
theta = data.theta(index);
Phi = data.phiSet{index}; 

M = data.M;



t = Phi{1};
phiAll = Phi{2};
if M.plotFall
    tFall = Phi{4};
    PhiFall = Phi{5};
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



xlabel("$q_1$",'interpreter','latex')
ylabel("$q_2$",'interpreter','latex')

end

