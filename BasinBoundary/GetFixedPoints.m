function [M] = GetFixedPoints(M)
%GETATTRACTORSET 

    %load data
    folder = "Data/FixedPoints";
if M.rhsString == "Duffing"

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
            data.Auto.FixedPoints{i}.W;
        end
        if foundFlag == false
            error("did not find fixedpoint data with that frequency")
        end
        M.Mrhs.FixedPoints =  FixedPointsClass(fpData.FP,fpData.phi,fpData.Stability,fpData.solution)
    else
        error("Parameters do not match with requested file");
    end
elseif M.rhsString == "TwoDuffing"
    if M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.Mrhs.kc == .1
        data = load(fullfile(folder,"FixedPointsTwoDuffing.mat"));
    elseif M.Mrhs.a1 == 1 & M.Mrhs.a3 == .3 & M.Mrhs.nu == .1 & M.Mrhs.F == .4 & M.Mrhs.kc == 0 & M.Mrhs.w == 1.4
        data1 = load(fullfile(folder,"FixedPoints.mat"));
        fpData1 = data1.Auto.FixedPoints{121}
        data2 =  load(fullfile(folder,"FixedPointsTwoDuffing.mat"));
        fpData2 = data2.Auto.FixedPoints{91};
        FP1 = fpData1.FP;
        FP2(1,:) = [FP1(1,:) FP1(1,:)]; %low amplitude mode 
        FP2(2,:) = [FP1(2,:) FP1(2,:)]; %saddle 1
        FP2(3,:) = [FP1(3,:) FP1(3,:)]; %high amplitude mode 
        FP2(4,:) = [FP1(3,:) FP1(1,:)]; %Localized 1 amplitude mode 
%         FP2(5,:) = [FP1(1,:) FP(1,:)]; %saddle 2 amplitude mode 
%         FP2(6,:) = [FP1(1,:) FP(1,:)]; %saddle 3 amplitude mode 
        FP2(5,:) = [FP1(1,:) FP1(3,:)]; %Localized 2 amplitude mode 
        fpData2.FP = FP2;
        fpData2.phi(6:7) =[];
        fpData2.stability(5) = 1;
        fpData2.stability(6:7) = [];   
        sol = fpData1.solution;
        t = sol{1}(:,1);
        for i = 1:length(sol)
            sol1{i} = [t,interp1(sol{i}(:,1),sol{i}(:,2),t),interp1(sol{i}(:,1),sol{i}(:,3),t)];
        end
        ii = [1 1];
        sol2{1} = [sol1{ii(1)} sol1{ii(2)}(:,2:3)];
        ii = [2 2];
        sol2{2} = [sol1{ii(1)} sol1{ii(2)}(:,2:3)];
        ii = [3 3];
        sol2{3} = [sol1{ii(1)} sol1{ii(2)}(:,2:3)];
        ii = [3 1];
        sol2{4} = [sol1{ii(1)} sol1{ii(2)}(:,2:3)];
        ii = [1 3];
        sol2{5} = [sol1{ii(1)} sol1{ii(2)}(:,2:3)];
        fpData2.solution = sol2;
        M.Mrhs.FixedPoints =  FixedPointsClass(fpData2.FP,fpData2.phi,fpData2.stability,fpData2.solution);
        return
    end

    %validate parameters
    m = data.Auto.M;
    if M.Mrhs.a1 == m.a & M.Mrhs.a3 == m.beta & M.Mrhs.nu == m.nu & M.Mrhs.F == m.F & M.Mrhs.kc == .1
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
        M.Mrhs.FixedPoints =  FixedPointsClass(fpData.FP,fpData.phi,fpData.stability,fpData.solution)
    else
        error("Parameters do not match with requested file");
    end
end

end


