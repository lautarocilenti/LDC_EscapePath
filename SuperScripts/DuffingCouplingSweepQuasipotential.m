function [] = DuffingCouplingSweepQuasipotential()
%FREQUENCYSWEEPQUASIPOTENTIAL Summary of this function goes here
cd ..
AddAllPaths()

Wlist = [1.4];
kcList = [0.01:0.01:0.06];
for i = 1:length(kcList)
    w = Wlist(1);
    kc = kcList(i);
     T = 2*pi/w;  dt = T/32; tf = 500*T; dT = T/2;
     tspan = [0:dt:tf];
    [data] = Main({'MRHS_w','MRHS_T','dt','tspan','tf','dT','MRHS_kc'},{w,T,dt,tspan,tf,dT,kc})
    Kc(i) = kc;
    Q(i) = data.minS;
    Data{i} = data;
    M = data.M;
    fileName = sprintf("Data/ActionPlot/CQ_%s_%s_kc%.2f.mat",M.rhsString,M.paramNote,kc);
    save(fileName,'data','kc','-v7.3');
end
M = data.M;
dateLog = datenum(datetime('now'));
formatOut = 'yyyy_mm_dd_HH_MM_SS_FFF';
fileName = sprintf("Data/ActionPlot/CQ_%s_%s_%s.mat",M.rhsString,M.paramNote,datestr(dateLog,formatOut));
save(fileName,'Data','Q','Kc','dateLog','-v7.3');
end

