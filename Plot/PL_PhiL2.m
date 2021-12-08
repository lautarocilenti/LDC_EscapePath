function [] = PL_PhiL2(phi,FP,T,M,fullpath)

if nargin ==4
    fullpath = false;
end

if fullpath
   phi = GenerateFullPath(phi,M); 
end



%PL_PHIL2 
jList = [1,4];
legendstr = {};
color = abs([rand() rand() rand()]-.5);
color = [ 0 0 0];

lspec = ["-",":"];
tqAll = [];
names = {"Escape Path","Trajectory"};
plotfall = true;
if plotfall
    jL = length(jList);
else
    jL = 1;
end
for j = 1:jL
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
    
    Path_L2 = vecnorm(x(:,1:M.dim),2,2);
    Path_L2q = vecnorm(xq,2,2);%phi index, 2nd index for path, all timesteps, first two columns
    plot(t,Path_L2,lspec(j),'Color',color,'DisplayName',sprintf('%s Continuous',names{j}),'linewidth',.5)
    hold on
    plot(tq,Path_L2q,'x','Color',color./2,'DisplayName',sprintf('%s Poincare',names{j}),'markersize',10,'linewidth',3)


    tOld = tq(end); xqOld = xq(end,:); tqAll = [tqAll tq];
end

for i = 1:size(FP,1)
    FP_L2 = norm(FP(i,:));
    if contains(M.attractorNames{i},"Saddle")
       marker = "p";
       color = [1 0 0];
       if M.dim > 6
           continue
       end
    else
        marker = "o";
        color = [0 .7 0];
    end
    s1 = scatter(tqAll,FP_L2*ones(size(tqAll)),10,marker,'DisplayName',sprintf('%s Poincare',M.attractorNames{i}),'markerfacecolor',color)
    s1.MarkerFaceAlpha = .7;
    s1.MarkerEdgeColor = color;
    hold on
end


xlim([0 max(tqAll)])
legend('location','southwest','fontsize',12)
xlabel('time (s)','fontsize',12)
ylabel('$||\textbf{q}||$','interpreter','latex','fontsize',12)
hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
end

