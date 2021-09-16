function [] = FrequencySweepQuasipotential()
%FREQUENCYSWEEPQUASIPOTENTIAL Summary of this function goes here


Wlist = [1.30:.01:1.55];
for i = 1:length(Wlist)
    w = Wlist(i);
     T = 2*pi/w;  dt = T/32; tf = 300*T;
     tspan = [0:dt:tf];
    [data] = Main({'MRHS_w','MRHS_T','dt','tspan','tf'},{w,T,dt,tspan,tf})
    W(i) = w;
    Q(i) = data.minS;
    Data{i} = data;
end
M = data.M;
dateLog = datenum(datetime('now'));
formatOut = 'yyyy_mm_dd_HH_MM_SS_FFF';
fileName = sprintf("Data/Frequency/FQ_%s_%s_%s.mat",M.rhsString,M.paramNote,datestr(dateLog,formatOut));
save(fileName,'Data','Q','W','dateLog','-v7.3');
end

