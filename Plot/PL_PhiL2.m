function [] = PL_PhiL2(phi,FP,T,M)
%PL_PHIL2 
jList = [1,4];
legendstr = {};
color = abs([rand() rand() rand()]-.5);
lspec = ["-",":"];
tqAll = [];
for j = 1:length(jList)
    xq = []; tq = [];
    t = phi{jList(j)};
    x = phi{jList(j)+1};
    tq = [0:T:max(t)];
    tq = tq(tq>=min(t));
    for d = 1:M.dim
        xq(:,d) = interp1(t,x(:,d),tq);
    end    

    if j == 2
       tq = [tOld tq];
       xq = [xqOld;xq];
    end
    
    Path_L2 = vecnorm(x,2,2);
    Path_L2q = vecnorm(xq,2,2);%phi index, 2nd index for path, all timesteps, first two columns
    plot(t,Path_L2,lspec(j),'Color',color,'DisplayName','Continuous')
    hold on
    plot(tq,Path_L2q,'x','Color',color,'DisplayName','Poincare')


    tOld = tq(end); xqOld = xq(end,:); tqAll = [tqAll tq];
end



for i = 1:size(FP,1)
    FP_L2 = norm(FP(i,:));
    plot(tqAll,FP_L2*ones(size(tqAll)),'.','DisplayName','FP Poincare','markersize',10)
    hold on
end

xlim([0 max(tqAll)])
end

