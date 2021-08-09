function [phiSetOut,msLog] = MinSearch1DLocalRegion(phiSet,msLog,nLM,M,regionCenters)
%MINSEARCH1D                 
% Identify all local minima (all nearest neighbors have larger energy)
% Select up to nLM local minima with smallest energy
% Generate initial conditions between the neighbors to the local minima and the local minima.
% Run rise and fall method for new initial conditions

theta = msLog{end}{1};
S = msLog{end}{2};

%Augmented S
thetaAug = [theta theta(1)];
sAug = [S;S(1)]; %adding the first term to the end to complete the circle



%Neighbor Differentials
ds1 = [0 ;diff(sAug)]; 
ds2 = [diff(sAug);0];

iLMSet = find(ds1<0 & ds2>0); %index of local minima
eps = .5;
counter = 1;
for iRegion = 1:length(regionCenters)
    rc = regionCenters(iRegion); 
    iThetaRegion = find(thetaAug < rc+eps & thetaAug > rc-eps);
    iLMRegion = find(ismember(iLMSet,iThetaRegion));
    [~,iSort] = sort(S(iLMRegion),'ascend');
    iLMRegion = iLMRegion(iSort);
    L = length(iLMRegion);
    if L > 3
        iLM(counter:counter+3) = iLMSet(iLMRegion(1:4));  
        counter = counter + 4;
    elseif L > 2
        iLM(counter:counter+2) = iLMSet(iLMRegion(1:3));  
        counter = counter + 3;
    elseif L > 1
        iLM(counter:counter+1) = iLMSet(iLMRegion(1:2));  
        counter = counter + 2;
    else
        iLM(counter) = iLMSet(iLMRegion(1));  
        counter = counter +1;
    end
end



sLM = S(iLM); %energy of local minima

%Limit to max number of local minima to use
if length(sLM) > nLM
    [~,iSort] = sort(sLM,'ascend');
    iLM = iLM(iSort);
    iLM = iLM(1:nLM);
    sLM = S(iLM);
end




%generate new initial conditions
iLeft = iLM - 1; iRight = iLM +1;

thetaNew =  mean([thetaAug(iLM) thetaAug(iLM); thetaAug(iLeft) thetaAug(iRight)]);

%Run Rise and Fall Method for new initial conditions
[phiSetNew,SNew] = GetNewPhiSets(thetaNew,M);

%store output data
thetaOut = [theta thetaNew];
phiSetOut = [phiSet phiSetNew];
Sout = [S ;SNew];

[thetaOut,iSortTheta] = sort(thetaOut,'ascend');
phiSetOut = phiSetOut(iSortTheta);
Sout = Sout(iSortTheta);

msLog{end+1} = {thetaOut,Sout};

end

function [phiSetNew,SNew] = GetNewPhiSets(theta,M)
     phiSetNew = PostProcessTrajectories2(IntegrateRHSinPeriods(GenerateInitialConditionsFloquet(theta,M),M),M);
     SNew = IntegrateLagrangian(phiSetNew,M);
end


