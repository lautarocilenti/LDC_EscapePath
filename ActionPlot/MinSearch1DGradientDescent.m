function [phiSetOut,msLog,TerminateFlag] = MinSearch1D(phiSet,msLog,nLM,M)
%MINSEARCH1D                 
% Identify all local minima (all nearest neighbors have larger energy)
% Select up to nLM local minima with smallest energy
% Generate initial conditions between the neighbors to the local minima and the local minima.
% Run rise and fall method for new initial conditions



theta = msLog{end}{1};
S = msLog{end}{2};
DescentPrev = msLog{end}{3};




if DescentPrev.Count == 0

    [thetaCurrent,sCurrent,~] = IdentifyLocalMinima(theta,S,nLM);
    
    runDescentOnTheta = true(size(thetaCurrent));
    descentStep = M.descent.Gamma*ones(size(thetaCurrent));
    sPrev = sCurrent;
    fdCostPrev = NaN(size(thetaCurrent));
    iRepeat = false(size(thetaCurrent));
    iDisc = false(size(thetaCurrent));

    
else
    c = DescentPrev.Count;
    thetaCurrent = DescentPrev.thetaNew(c,:);
    sCurrent =  DescentPrev.sNew(c,:);
    runDescentOnTheta =  DescentPrev.runDescentOnTheta(c,:);
    sPrev =  DescentPrev.sCurrent(c,:);
    fdCostPrev = DescentPrev.fdCost(c,:);
    descentStep = DescentPrev.descentStep(c,:);
    iRepeat = DescentPrev.iRepeatNext;
    iDisc = DescentPrev.iDisc;
    

    
end



[fdCost] = FiniteDifferenceCentralGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M)

iStop = find(abs(fdCost)<1E-2 | abs(fdCost)>1E100); %terminate continuation of descent
runDescentOnTheta(iStop) = 0;
% 


thetaNew = thetaCurrent;
igradientTheta = runDescentOnTheta & ~iDisc;
thetaNew(igradientTheta) =  thetaCurrent(igradientTheta)-descentStep(igradientTheta).*fdCost(igradientTheta); %calculate new descent cost
iDiscTheta = runDescentOnTheta & iDisc;
thetaNew(iDiscTheta) =  thetaCurrent(iDiscTheta)-descentStep(iDiscTheta).*sign(fdCost(iDiscTheta)); %calculate new descent cost





if DescentPrev.Count > 0
    [thetaNew,descentStep,runDescentOnTheta,iDisc] = AdaptiveStepCheck(theta,S,thetaCurrent,thetaNew,fdCost,descentStep,runDescentOnTheta,sPrev,iDisc,M);
end

figure(100)
xx = linspace(0,2*pi,100);
yy = interp1(theta,S,xx,'linear');
plot(xx,yy,'.-')
hold on

plot(thetaCurrent,sCurrent,'o')
plot(thetaNew,sCurrent,'x')
drawnow()
hold off

if any(runDescentOnTheta)
    %Run Rise and Fall Method for new initial conditions
    [sNew,phiSetNew] = CostFunction(thetaNew(runDescentOnTheta),M);
    
    %store output data
    thetaOut = [theta thetaNew(runDescentOnTheta)];
    phiSetOut = [phiSet phiSetNew];
    sOut = [S  sNew];
    sNewOut = sPrev; sNewOut(runDescentOnTheta) = sNew;

    [thetaOut,iSortTheta] = sort(thetaOut,'ascend');
    phiSetOut = phiSetOut(iSortTheta);
    sOut = sOut(iSortTheta);
    TerminateFlag = false;
else
    thetaOut = theta;
    phiSetOut = phiSet;
    sOut = S;
    sNewOut = sPrev;
    TerminateFlag = true;
end

%Allow cancellation of bad moves
iCancelLastMove = (sNewOut>=1.2*sPrev); %indices that increased the cost by more than 20% 
thetaNew(iCancelLastMove) = thetaCurrent(iCancelLastMove);
sNewOut(iCancelLastMove) = sPrev(iCancelLastMove);
fdCost(iCancelLastMove) = fdCostPrev(iCancelLastMove);
iRepeatNext = false(size(thetaCurrent));
iRepeatNext(iCancelLastMove) = true;

%Allow step change given not too bad moves
iStepChange = (sNewOut>sPrev); %indices that increased the cost 
iStepChange = iStepChange & ~iCancelLastMove;
descentStep(iStepChange) = descentStep(iStepChange)/2;


Descent = DescentPrev;
Descent.Count = DescentPrev.Count+1;
c = Descent.Count;

Descent.runDescentOnTheta(c,:) = runDescentOnTheta;
Descent.theta{c} = thetaOut;
Descent.S{c} = sOut;
Descent.thetaCurrent(c,:) = thetaCurrent;
Descent.sCurrent(c,:) = sCurrent;
Descent.thetaNew(c,:) = thetaNew;
Descent.sNew (c,:) = sNewOut;
Descent.fdCost(c,:) = fdCost;
Descent.descentStep(c,:) = descentStep;
Descent.iRepeatNext = iRepeatNext;
Descent.iDisc = iDisc;
    

msLog{end+1} = {thetaOut,sOut,Descent};


end

function [thetaLM,sLM,iLM] = IdentifyLocalMinima(theta,S,nLM)

%Augmented optimization space
    thetaAug = [theta theta(1)];
    sAug = [S S(1)]; %adding the first term to the end to complete the circle

 %Find local minima
    %Neighbor Differentials
    ds1 = [0 diff(sAug)]; 
    ds2 = [diff(sAug) 0];


    iLM = find(ds1<0 & ds2>0); %index of local minima
    sLM = S(iLM); %energy of local minima

    %Limit to max number of local minima to use
    if length(sLM) > nLM
        [~,iSort] = sort(sLM,'ascend');
        iLM = iLM(iSort);
        iLM = iLM(1:nLM);
        sLM = S(iLM);
    end
    thetaLM = thetaAug(iLM);
end

function [fdCost] = FiniteDifferenceGradientRightOnly(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,M)
    fdStep = M.descent.fdStep;
    iFindFDCost = find(~iRepeat & runDescentOnTheta & ~iDisc); %check which indices need new gradient
    fdCost = fdCostPrev; %load prior gradient values
    if any(iFindFDCost)
        fdTheta = [thetaCurrent(iFindFDCost)+fdStep]; %finite difference calculation angles
        [fdS,~] = CostFunction(fdTheta,M); %cost near initial point 
        fdSInit = [sCurrent(iFindFDCost)]; %initial cost
        fdGradient = (fdS - fdSInit)./fdStep; %slope of change in cost
        fdCost(iFindFDCost)  = fdGradient; %finite difference in cost 
    end
end
% 
function [fdCost] = FiniteDifferenceCentralGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M)
    fdStep = M.descent.fdStep*ones(size(thetaCurrent));
    fdStep(iDisc) = descentStep(iDisc);
    iFindFDCost = find(~iRepeat & runDescentOnTheta); %check which indices need new gradient
    fdCost = fdCostPrev; %load prior gradient values
    if any(iFindFDCost)
        fdStepVector = [fdStep(iFindFDCost) -fdStep(iFindFDCost)];
        fdTheta = [thetaCurrent(iFindFDCost) thetaCurrent(iFindFDCost)]+fdStepVector ; %finite difference calculation angles
        [fdS,~] = CostFunction(fdTheta,M); %cost near initial point 
        fdSInit = [sCurrent(iFindFDCost) sCurrent(iFindFDCost)]; %initial cost
        fdGradient = (fdS - fdSInit)./fdStepVector; %slope of change in cost
        fdGradient = [fdGradient(1:length(fdGradient)/2);fdGradient(length(fdGradient)/2+1:end)]; %reshape array 
        fdCentralGradient = sum(fdGradient,1)/2; % calculate the central gradient (mean of the two gradients)
        fdCost(iFindFDCost)  = fdCentralGradient; %finite difference in cost 
    end
end

function [thetaNew,descentStep,runDescentOnTheta,iDisc] = AdaptiveStepCheck(theta,S,thetaCurrent,thetaNew,fdCost,descentStep,runDescentOnTheta,sPrev,iDisc,M)

    iCheck = SmartChecks1Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sPrev);
    iDisc = iDisc | iCheck;
    %adaptive step to account for smart checks
    if any(iCheck)

       while(any(iCheck))
           iNeedsAdapt = find(iCheck);
           for j = 1:length(iNeedsAdapt)
               jAdapt = iNeedsAdapt(j);
               if descentStep(jAdapt) > M.descent.minGamma %limit to the adaptive step, if it shrinks too much end this descent
                 descentStep(jAdapt) = descentStep(jAdapt)./1.01;
                 thetaNew(jAdapt) =  thetaCurrent(jAdapt)-descentStep(jAdapt).*sign(fdCost(jAdapt));
               else
                   fprintf("Terminating due to small adaptive step\n")
                   thetaNew(jAdapt) = thetaCurrent(jAdapt);
                   runDescentOnTheta(jAdapt) = 0;
               end
           end
           iCheck = SmartChecks1Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sPrev);
       end
%        fprintf("Using an adaptive step, smallest step = %e\n", min(descentStep))
    end
    if any(ismember(thetaNew(runDescentOnTheta),theta))
            error("New theta is the same as old theta\n")
    end
end


function [SNew,phiSetNew] = CostFunction(theta,M)
    M.progressbar = false;
     phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
     SNew = IntegrateLagrangian(phiSetNew,M); 
end

function [iCheck] = SmartChecks1Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sOld)

        for i = 1:length(thetaNew)
            SPredict(1,i) = interp1(theta,S,thetaNew(i),'nearest');
        end
        iCheck2 = ismember(thetaNew,theta);

        iCheck = SPredict>=1.2*sOld;
        iCheck = (iCheck | iCheck2) & runDescentOnTheta ;


end

