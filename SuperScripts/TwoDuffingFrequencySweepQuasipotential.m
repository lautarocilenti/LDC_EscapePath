function [] = FrequencySweepQuasipotential()
%FREQUENCYSWEEPQUASIPOTENTIAL Summary of this function goes here
cd ..
AddAllPaths()
data = load("FixedPointsDuffing.mat");
Wlist = cellfun(@(x) x.W,data.Auto.FixedPoints);
Wlist(1) = [];
Wlist = [1.35:.01:1.54];
for i = 1:2%length(Wlist)
    w = Wlist(i);
    if w<1.4
        rA1 = 1;
    else
        rA1 = .1;
    end
     T = 2*pi/w;  dt = T/32; tf = 500*T; dT = T/2;
     tspan = [0:dt:tf];
    [data] = Main({'MRHS_w','MRHS_T','dt','tspan','tf','dT','rA1'},{w,T,dt,tspan,tf,dT,rA1})
    W(i) = w;
    Q(i) = data.minS;
    Data{i} = data;
    M = data.M;
    fileName = sprintf("Data/ActionPlot/FQ_%s_%s_w%f.mat",M.rhsString,M.paramNote,w);
    save(fileName,'data','w','-v7.3');
end
M = data.M;
dateLog = datenum(datetime('now'));
formatOut = 'yyyy_mm_dd_HH_MM_SS_FFF';
fileName = sprintf("Data/ActionPlot/FQ_%s_%s_%s.mat",M.rhsString,M.paramNote,datestr(dateLog,formatOut));
save(fileName,'Data','Q','W','dateLog','-v7.3');
end

