function [phiSetOut,msLog,TerminateFlag] = MinSearch1D(phiSet,msLog,nLM,M)
%MINSEARCH1D                 
% Identify all local minima (all nearest neighbors have larger energy)
% Select up to nLM local minima with smallest energy
% Generate initial conditions between the neighbors to the local minima and the local minima.
% Run rise and fall method for new initial conditions

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
    
    thetai = thetaAug(iLM);

    aCounter = ones(size(thetai));
    
else
    thetai = DescentPrev.thetaNew;
    sStart = DescentPrev.SNew;
    thetaOld = DescentPrev.thetaStart;
    sOld = DescentPrev.sStart; 
    sOld(DescentPrev.iStop) = [];
    aCounter = DescentPrev.aCounter;
    iCancelLastMove = find(sStart>=1.2*sOld); %indices that increased the cost by more than 20% 
    thetai(iCancelLastMove) = thetaOld(iCancelLastMove);
    sStart(iCancelLastMove) = sOld(iCancelLastMove);
   
    
end


fd_step = .01; %finite difference step size
C_gamma = .5; %gradient descent initial step size

Sright = CostFunction(thetai+fd_step,M); %cost near initial point
fdCost  = (Sright-sStart)/fd_step; %finite difference in cost

iStop = find(abs(fdCost)<5E-2 | abs(fdCost)>1E1); %terminate continuation of descent

if length(iStop) >0 
   fprintf("Terminating descent for %d paths \n" ,length(iStop)) 
end

FD_Copy = fdCost %copy fd cost
thetai(iStop) = [];  fdCost(iStop) = []; %clear terminated descent paths
aCounter(iStop) = [];
gamma = C_gamma./(2*aCounter);
thetaNew =  thetai-gamma.*fdCost'; %calculate new descent cost

iRepeat = ismember(thetaNew,thetaAug);
%adaptive step to avoid repeats
if any(iRepeat)
   fprintf("Using an adaptive step, counter = %d\n", max(aCounter))
   while(any(iRepeat))
       ii = find(iRepeat);
       thetaStart = thetai(iRepeat);
       fdStart(1,:) = fdCost(iRepeat);
       thetaNew(iRepeat) = []; %clear repeats
       for j = 1:length(ii)
           if aCounter(ii(j)) < 10 %limit to the adaptive step, if it shrinks too much end this descent
               for i = 1:length(thetaStart)
                 thetaNew(end+1) =  thetaStart(i)-C_gamma./(2*aCounter(ii(j)))*fdStart(i);
               end
           else
               fprintf("Terminating due to small adaptive step\n")
               iStop = [iStop ii(j)];
           end
       end
       iRepeat = ismember(thetaNew,thetaAug);
       aCounter(ii(j)) = aCounter(ii(j))+1;
   end
end



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
Descent.thetaStart = thetai;
Descent.sStart = sStart;
Descent.Count = DescentPrev.Count+1;
Descent.iStop = iStop;
Descent.aCounter = aCounter;


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

