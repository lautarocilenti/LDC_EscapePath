function [phiSetOut,msLog,TerminateFlag] = MinSearch1D(phiSet,msLog,nLM,M)
%MINSEARCH1D                 
% Identify all local minima (all nearest neighbors have larger energy)
% Select up to nLM local minima with smallest energy
% Generate initial conditions between the neighbors to the local minima and the local minima.
% Run rise and fall method for new initial conditions


%clean up this code specially the termination of paths

theta = msLog{end}{1};
S = msLog{end}{2};
DescentPrev = msLog{end}{3};

%Augmented S
thetaAug = [theta theta(1)];
sAug = [S;S(1)]; %adding the first term to the end to complete the circle

if DescentPrev.Count == 0
    %Find local minima


    %Neighbor Differentials
    ds1 = [0 ;diff(sAug)]; 
    ds2 = [diff(sAug);0];


    iLM = find(ds1<0 & ds2>0); %index of local minima
    sStart = S(iLM); %energy of local minima

    %Limit to max number of local minima to use
    if length(sStart) > nLM
        [~,iSort] = sort(sStart,'ascend');
        iLM = iLM(iSort);
        iLM = iLM(1:nLM);
        sStart = S(iLM);
    end
    
    thetaInitial = thetaAug(iLM);
    fdCostPrev = NaN(size(thetaInitial));
    iTheta = 1:length(thetaInitial);

    aCounter = ones(size(thetaInitial));
    
else
    %load current step
    thetaInitial = DescentPrev.thetaNew;
    sStart = DescentPrev.SNew;
    iTheta = DescentPrev.iTheta;
    %load previous step
    thetaOld = DescentPrev.thetaStart;
    sOld = DescentPrev.sStart; 
    sOld(DescentPrev.iStop) = [];
    aCounter = DescentPrev.aCounter;
    %determine if cancelling the previous step
    iCancelLastMove = find(sStart>=1.2*sOld); %indices that increased the cost by more than 20% 
    thetaInitial(iCancelLastMove) = thetaOld(iCancelLastMove);
    sStart(iCancelLastMove) = sOld(iCancelLastMove);
    
    fdCostPrev = NaN(size(thetaInitial));
    fdCostPrev(iCancelLastMove) = DescentPrev.fdCost(iCancelLastMove);
   
    
end

C_gamma = 1; %gradient descent initial step size
for i = 1:length(aCounter)
    fd_step(i,1) = .01;
end


iFindFDCost = isnan(fdCostPrev);
fdCost = fdCostPrev;
if any(iFindFDCost)
    
    Sright = CostFunction(thetaInitial(iFindFDCost)+fd_step(iFindFDCost)',M); %cost near initial point
    fdCost(iFindFDCost)  = (Sright-sStart(iFindFDCost))./fd_step(iFindFDCost); %finite difference in cost
end

iStop = find(abs(fdCost)<1E-2 | abs(fdCost)>1E1); %terminate continuation of descent

if length(iStop) >0 
   fprintf("Terminating descent for %d paths \n" ,length(iStop)) 
end

FD_Copy = fdCost %copy fd cost

thetaInitial(iStop) = [];  fdCost(iStop) = []; iTheta(iStop) = []; %clear terminated descent paths
aCounter(iStop) = [];  
gamma = C_gamma./(2*aCounter);
thetaNew =  thetaInitial-gamma.*fdCost; %calculate new descent cost

if DescentPrev.Count > 0
    sOld_Copy = sOld; sOld_Copy(iStop) = [];
% iRepeat = ismember(thetaNew,thetaAug);
    iCheck = SmartChecks(theta,S,thetaNew,sOld_Copy);
    %adaptive step to account for smart checks
    if any(iCheck)

       while(any(iCheck))
           ii = find(iCheck);
           thetaStart = thetaInitial(iCheck);
           fdStart(1,:) = fdCost(iCheck);
           thetaNew(iCheck) = []; %clear repeats
           for j = 1:length(ii)
               if aCounter(ii(j)) < 100000 %limit to the adaptive step, if it shrinks too much end this descent
                 thetaNew(1,end+1) =  thetaStart(j)-C_gamma./(2.*aCounter(ii(j)))*fdStart(j);
               else
                   fprintf("Terminating due to small adaptive step\n")
                   iStop = [iStop ii(j)];
               end
           end
           iCheck = SmartChecks(theta,S,thetaNew,sOld_Copy);
           aCounter(ii) = aCounter(ii)+1;
       end
       fprintf("Using an adaptive step, counter = %d\n", max(aCounter))
    end
end
%Check closest point that 



if ~isempty(thetaNew)
    %Run Rise and Fall Method for new initial conditions
    [phiSetNew,SNew] = GetNewPhiSets(thetaNew,M);

    %store output data
    thetaOut = [theta thetaNew];
    phiSetOut = [phiSet phiSetNew];
    Sout = [S ; SNew];

    [thetaOut,iSortTheta] = sort(thetaOut,'ascend');
    phiSetOut = phiSetOut(iSortTheta);
    Sout = Sout(iSortTheta);
    TerminateFlag = false;
else
    thetaOut = theta;
    phiSetOut = phiSet;
    Sout = S;
    SNew = [];
    TerminateFlag = true;
end


Descent.ithetaNew = find(ismember(thetaOut,thetaNew));
Descent.thetaNew = thetaNew;
Descent.SNew = SNew;
Descent.thetaStart = thetaInitial;
Descent.sStart = sStart;
Descent.Count = DescentPrev.Count+1;
Descent.iStop = iStop;
Descent.aCounter = aCounter;
Descent.fdCost = fdCost;
Descent.iTheta = iTheta;

msLog{end+1} = {thetaOut,Sout,Descent};


end



function [phiSetNew,SNew] = GetNewPhiSets(theta,M)
     phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
     SNew = IntegrateLagrangian(phiSetNew,M); 
end

function SNew = CostFunction(theta,M)
     phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
     SNew = IntegrateLagrangian(phiSetNew,M); 
end

function [iCheck] = SmartChecks(theta,S,thetaNew,sOld)
    if isempty(thetaNew)
        iCheck = [];
    else
        for i = 1:length(thetaNew)
            SPredict = interp1(theta,S,thetaNew(i),'nearest');
        end
        iCheck = SPredict>=1.2*sOld;
        iCheck = iCheck';
    end


end

