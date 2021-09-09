function [] = PL_L2MPEPPath(data)
%PL_L2PATH Summary of this function goes here
%   Detailed explanation goes here
M = data.M;



theta = data.theta(data.minPhiIndex);
minPhi = data.phiSet{data.minPhiIndex}; 
minS = data.minS;

FP = M.Mrhs.FixedPoints.FP;
t = M.tspan;
T = 2*pi/M.Mrhs.w;
tmax = 0;

phi = minPhi;

jList = [1];
legendstr = {};
color = [rand() rand() rand()];
lspec = ["-","--"];
for j = 1:length(jList)
    t = phi{jList(j)};
    x1 = phi{jList(j)+1}(:,1);
    x2 = phi{jList(j)+1}(:,2);
    x3 = phi{jList(j)+1}(:,3);
    x4 = phi{jList(j)+1}(:,4);
    tq = [min(t):T:max(t)];
    x1q = interp1(t,x1,tq)';
    x2q = interp1(t,x2,tq)';
    x3q = interp1(t,x3,tq)';
    x4q = interp1(t,x4,tq)';
    if j == 2
       tq = [tOld tq];
       x1q=[x1Old;x1q]; x2q=[x2Old;x2q];
    end
    Path_L2 = vecnorm([x1,x2,x3,x4],2,2);%phi index, 2nd index for path, all timesteps, first two columns
    Path_L2q = vecnorm([x1q,x2q,x3q,x4q],2,2);%phi index, 2nd index for path, all timesteps, first two columns
    plot(t,Path_L2,lspec(j),'Color',color,'DisplayName','Continuous')
    hold on
    plot(tq,Path_L2q,'x','Color',color,'DisplayName','Poincare')

    if max(t) > tmax
       tmax = max(t); 
    end
    tOld = tq(end); x1Old = x1q(end); x2Old = x2q(end);
end



for i = 1:size(FP,1)
    FP_L2 = norm(FP(i,:));
    plot(tq,FP_L2*ones(size(tq)),'.','DisplayName','FP Poincare','markersize',10)
    hold on
end



end

