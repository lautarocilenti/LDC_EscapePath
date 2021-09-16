function [] = PL_L2Path(data)
%PL_L2PATH Summary of this function goes here
%   Detailed explanation goes here
M = data.M;
FP = M.Mrhs.FixedPoints.FP;
T = 2*pi/M.Mrhs.w;
tmax = 0;
phiSet = data.phiSet;
for i = 1:length(phiSet)
    if M.plotFall
            jList = [1,4];
    else
        jList = [1];
    end
    
    color = [rand() rand() rand()];
    lspec = ["-","--"];
    for j = 1:length(jList)
        t = phiSet{i}{jList(j)};
        x1 = phiSet{i}{jList(j)+1}(:,1);
        x2 = phiSet{i}{jList(j)+1}(:,2);
        x3 = phiSet{i}{jList(j)+1}(:,3);
        x4 = phiSet{i}{jList(j)+1}(:,4);
        tq = min(t):T:max(t);
        x1q = interp1(t,x1,tq)';
        x2q = interp1(t,x2,tq)';
        x3q = interp1(t,x3,tq)';
        x4q = interp1(t,x4,tq)';
        
        
        if j == 2
           tq = [tOld tq];
           x1q=[x1Old;x1q]; x2q=[x2Old;x2q]; x3q=[x3Old;x3q]; x4q=[x4Old;x4q];
        end
        Path_L2 = vecnorm([x1q,x2q,x3q,x4q],2,2);%phi index, 2nd index for path, all timesteps, first two columns
        plot(tq,Path_L2,lspec(j),'Color',color);
        hold on
        if max(t) > tmax
           tmax = max(t); 
        end
        tOld = tq(end);
        x1Old = x1q(end); x2Old = x2q(end); x3Old = x3q(end); x4Old = x4q(end); 
    end
end


for i = 1:size(FP,1)
    FP_L2 = norm(FP(i,:));
    plot([0,tmax],[FP_L2,FP_L2])
    hold on
end



end

