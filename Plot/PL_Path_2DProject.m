function [] = PL_Path_2DProject(data,index)
%PL_MPEP 
theta = data.theta(index);
Phi = data.phiSet{index}; 

M = data.M;




[Phi] = GenerateFullPath(Phi,M)

t = Phi{1};
phiAll = Phi{2};
if M.plotFall
    tFall = Phi{4};
    PhiFall = Phi{5};
end

hold on
for i = 1:M.dim/2
    subplot(1,M.dim/2,i)
    phi = phiAll(:,(i-1)*2+1:(i)*2);
        msize = 4;
        color = [rand() rand() rand()]*.5;
        [tq,xq] = InterpolateToPhaseAngle(t,phi,M);
        plot(xq(:,1),xq(:,2),'.-','markersize',msize,'Color',color,'displayname','Rise');
        hold on
        if M.plotFall
            phiFall = PhiFall(:,(i-1)*2+1:(i)*2);
            [tqFall,xqFall] = InterpolateToPhaseAngle(tFall,phiFall,M);
            plot([xq(end,1);xqFall(:,1)],[xq(end,2);xqFall(:,2)],'.--','markersize',msize,'Color',color,'displayname','Fall');      
        end
    PL_Attractors(data.attractors,(i-1)*2,M.Mrhs.FixedPoints.names);
    axis([-3  3 -3 3])
    xL = sprintf('q_%d',(i-1)*2+1);
    xlabel(['$',xL,'$'],'interpreter','latex')
    yL = sprintf('q_%d',(i-1)*2+2);
    ylabel(['$',yL,'$'],'interpreter','latex')
    legend("location","southoutside")
end




end

