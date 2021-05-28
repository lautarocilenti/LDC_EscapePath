function [] = PL_L2Path(data)
%PL_L2PATH Summary of this function goes here
%   Detailed explanation goes here
M = data.M;
FP = M.Mrhs.FixedPoints.FP;
T = 2*pi/M.Mrhs.w;
tmax = 0;
phiSet = data.phiSet;
for i = 1:length(phiSet)
    if strcmp(M.terminateType,"Saddle")
            jList = [1];
    else
        jList = [1,4];
    end
    
    color = [rand() rand() rand()];
    lspec = ["-","--"];
    for j = 1:length(jList)
        t = phiSet{i}{jList(j)};
        x1 = phiSet{i}{jList(j)+1}(:,1);
        x2 = phiSet{i}{jList(j)+1}(:,2);
        tq = min(t):T:max(t);
        x1q = interp1(t,x1,tq)';
        x2q = interp1(t,x2,tq)';
        if j == 2
           tq = [tOld tq];
           x1q=[x1Old;x1q]; x2q=[x2Old;x2q];
        end
        Path_L2 = vecnorm([x1q,x2q],2,2);%phi index, 2nd index for path, all timesteps, first two columns
        plot(tq,Path_L2,lspec(j),'Color',color);
        hold on
        if max(t) > tmax
           tmax = max(t); 
        end
        tOld = tq(end);
        x1Old = x1q(end); x2Old = x2q(end);
    end
end


for i = 1:size(FP,1)
    FP_L2 = norm(FP(i,:));
    plot([0,tmax],[FP_L2,FP_L2])
    hold on
end



end

