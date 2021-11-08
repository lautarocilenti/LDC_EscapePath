function [] = PL_Path_2DProject(data,index)
%PL_MPEP 
theta = data.theta(index);
Phi = data.phiSet{index}; 

M = data.M;




[Phi] = GenerateFullPath(Phi,M);

t = Phi{1};
phiAll = Phi{2};
if M.plotFall
    tFall = Phi{4};
    PhiFall = Phi{5};
end

hold on
color = [rand() rand() rand()]*.5;
for i = 1:M.dim/2
    subplot(M.dim/2,1,i)
    phi = phiAll(:,(i-1)*2+1:(i)*2);
        msize = 10;
        
        [tq,xq] = InterpolateToPhaseAngle(t,phi,M);
        plot(xq(:,1),xq(:,2),'.-','markersize',msize,'Color',color,'displayname','Escape Path');
        hold on
        if M.plotFall
            phiFall = PhiFall(:,(i-1)*2+1:(i)*2);
            [tqFall,xqFall] = InterpolateToPhaseAngle(tFall,phiFall,M);
            plot([xq(end,1);xqFall(:,1)],[xq(end,2);xqFall(:,2)],'.-.','markersize',msize,'Color',color,'displayname','Duffing Trajectory');      
        end
    if isfield(M.Mrhs,"FixedPoints")
        PL_Attractors(data.attractors,M.Mrhs.FixedPoints.names,(i-1)*2);
    else
        PL_Attractors(data.attractors,[],(i-1)*2);
    end
    axis([-2  2 -2 3])
    xL = sprintf('q_%d',(i-1)*2+1);
    xlabel(['$',xL,'$'],'interpreter','latex')
    yL = sprintf('q_%d',(i-1)*2+2);
    ylabel(['$',yL,'$'],'interpreter','latex')
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
    legend("location","eastoutside")
end




end

