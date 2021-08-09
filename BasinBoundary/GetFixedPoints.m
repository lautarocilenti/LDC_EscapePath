function [M] = GetFixedPoints(M)
%GETATTRACTORSET 


%load data
folder = "Data\FixedPoints";
if M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 
    data = load(fullfile(folder,"FixedPointsDuffing.mat"));
end

%validate parameters
m = data.Auto.M;
if M.Mrhs.a1 == m.a & M.Mrhs.a3 == m.beta & M.Mrhs.nu == m.nu & M.Mrhs.F == m.F
    foundFlag = false;
    for i = 1:length(data.Auto.FixedPoints)
        if round(data.Auto.FixedPoints{i}.W,8) == round(M.Mrhs.w,8)
            fpData = data.Auto.FixedPoints{i};
            foundFlag = true;
            break;
        end
        data.Auto.FixedPoints{i}.W
    end
    if foundFlag == false
        error("did not find fixedpoint data with that frequency")
    end
    M.Mrhs.FixedPoints =  FixedPointsClass(fpData.FP,fpData.phi,fpData.Stability,fpData.solution)
else
    error("Parameters do not match with requested file");
end


end


