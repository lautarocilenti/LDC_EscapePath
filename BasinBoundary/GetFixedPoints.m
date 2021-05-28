function [M] = GetFixedPoints(M)
%GETATTRACTORSET 


%load data
folder = "Data\FixedPoints";
if M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.Mrhs.w == 1.4
    data = load(fullfile(folder,"FixedPointsDuffing.mat"));
end

%validate parameters
m = data.Auto.M;
if M.Mrhs.a1 == m.a & M.Mrhs.a3 == m.beta & M.Mrhs.nu == m.nu & M.Mrhs.F == m.F & M.Mrhs.w == m.w
    fpData = data.Auto.FixedPoints;
    M.Mrhs.FixedPoints =  FixedPointsClass(fpData.FP,fpData.phi,fpData.Stability,fpData.solution)
else
    error("Parameters do not match with requested file");
end


end


