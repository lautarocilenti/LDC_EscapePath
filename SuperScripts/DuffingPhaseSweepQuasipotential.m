function [] = DuffingCouplingSweepQuasipotential()
%FREQUENCYSWEEPQUASIPOTENTIAL Summary of this function goes here
cd ..
AddAllPaths()

Wlist = [1.4];
ppList = [0:pi:pi];
for i = 1:length(ppList)
    w = Wlist(1);
    pp = ppList(i);
     T = 2*pi/w;  dt = T/32; tf = 500*T; dT = T/2;
     tspan = [0:dt:tf];
    [data] = Main({'MRHS_w','MRHS_T','dt','tspan','tf','dT','pp'},{w,T,dt,tspan,tf,dT,pp})
    PP(i) = pp;
    Q(i) = data.minS;
    Data{i} = data;
    M = data.M;
    fileName = sprintf("Data/ActionPlot/PQ_%s_%s_kc%.2f.mat",M.rhsString,M.paramNote,pp);
    save(fileName,'data','kc','-v7.3');
end
M = data.M;
dateLog = datenum(datetime('now'));
formatOut = 'yyyy_mm_dd_HH_MM_SS_FFF';
fileName = sprintf("Data/ActionPlot/PQ_%s_%s_%s.mat",M.rhsString,M.paramNote,datestr(dateLog,formatOut));
save(fileName,'Data','Q','Kc','dateLog','-v7.3');
end

