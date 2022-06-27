function [M] = GetFixedPoints(M)
%GETATTRACTORSET 

    %load data
    folder = "Data/FixedPoints";
    if M.methodTest
        M.attractors = [];
        M.attractorNames = [];
        M.Mrhs.iA = [];
        M.Mrhs.fA = [];
        return
    elseif M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.dim == 2
        
        data = load(fullfile(folder,"FixedPointsDuffing.mat"))
        m = data.Auto.M;
    elseif M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.Mrhs.kc == .1 & M.dim == 4
        data = load(fullfile(folder,"FixedPointsTwoDuffing.mat"));
    elseif M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.Mrhs.kc == .01 & M.dim == 4
         data = load(fullfile(folder,"FixedPointsTwoDuffingkcp01.mat"));
         m = data.Auto.M;
    elseif M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.dim == 6
        filename = sprintf("FixedPointsThreeDuffingkc%.2f.mat",M.Mrhs.kc);
        data = load(fullfile(folder,filename));
        m = data.Auto.FixedPoints{1}.M;
    elseif M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.dim == 10
        filename = sprintf("FixedPointsFiveDuffingkc%.2f.mat",M.Mrhs.kc);
        data = load(fullfile(folder,filename));
        m = data.Auto.FixedPoints{1}.M;
    elseif M.dim == 2
        filename = sprintf("OneCantilever_Coco.mat");
        data = load(fullfile(folder,filename));
        m = data.coco.M.Mrhs;
        coco = data.coco.Branchq;
        M.Mrhs.FixedPoints = CocoFixedPointsClass(coco,M.iAString,M.fAString,M.Mrhs.w);
        M.attractors = M.Mrhs.FixedPoints.FP;
        M.attractorNames = M.Mrhs.FixedPoints.names;
        M.Mrhs.iA = M.Mrhs.FixedPoints.iA;
        M.Mrhs.fA = M.Mrhs.FixedPoints.fA;
        return
    elseif M.dim == 4
        filename = sprintf("TwoCantileverCoco_kc0.002000.mat");
        data = load(fullfile(folder,filename));
        m = data.coco.M.Mrhs;
        coco = data.coco.Branchq;
        M.Mrhs.FixedPoints = CocoFixedPointsClass(coco,M.iAString,M.fAString,M.Mrhs.w);
        M.attractors = M.Mrhs.FixedPoints.FP;
        M.attractorNames = M.Mrhs.FixedPoints.names;
        M.Mrhs.iA = M.Mrhs.FixedPoints.iA;
        M.Mrhs.fA = M.Mrhs.FixedPoints.fA;
        return
    end

    %validate parameters
    
    if M.Mrhs.a1 == m.a & M.Mrhs.a3 == m.beta & M.Mrhs.nu == m.nu & M.Mrhs.F == m.F
        foundFlag = false;
        for i = 1:length(data.Auto.FixedPoints)
            if round(data.Auto.FixedPoints{i}.W,8) == round(M.Mrhs.w,8)
                fpData = data.Auto.FixedPoints{i};
                foundFlag = true;
                break;
            end
            data.Auto.FixedPoints{i}.W;
        end
        if foundFlag == false
            error("did not find fixedpoint data with that frequency")
        end
        if isfield(fpData,'L2')
            M.Mrhs.FixedPoints =  FixedPointsClass(fpData.FP,fpData.phi,fpData.stability,fpData.solution,M.iAString,M.fAString,fpData.L2);
        else
             M.Mrhs.FixedPoints =  FixedPointsClass(fpData.FP,fpData.phi,fpData.stability,fpData.solution,M.iAString,M.fAString);
        end
        M.attractors = M.Mrhs.FixedPoints.FP;
        M.attractorNames = M.Mrhs.FixedPoints.names;
        M.Mrhs.iA = M.Mrhs.FixedPoints.iA;
        M.Mrhs.fA = M.Mrhs.FixedPoints.fA;
    else
        error("Parameters do not match with requested file");
    end

end


