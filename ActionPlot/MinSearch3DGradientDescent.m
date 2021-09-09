function [phiSetOut,msLog,TerminateFlag] = MinSearch3DGradientDescent(phiSet,msLog,nLM,M)
%MINSEARCH1D                 
% Identify all local minima (all nearest neighbors have larger energy)
% Select up to nLM local minima with smallest energy
% Generate initial conditions between the neighbors to the local minima and the local minima.
% Run rise and fall method for new initial conditions


%extract data from min search log
theta = msLog{end}{1};
S = msLog{end}{2};
DescentPrev = msLog{end}{3};




if DescentPrev.newStart  %first loop

    [~,sCurrent,iLM] = IdentifyLocalMinima(vecnorm(theta,2,1),S,nLM);
    thetaCurrent = theta(:,iLM);
    runDescentOnTheta = true(size(sCurrent));
    descentStep = M.descent.Gamma*ones(size(thetaCurrent));
    sPrev = sCurrent;
    fdCostPrev = NaN(size(thetaCurrent));
    iRepeat = false(size(sCurrent));
    iDisc = false(size(sCurrent));

    
else
    c = DescentPrev.Count;
    thetaCurrent = DescentPrev.thetaNew{c};
    sCurrent =  DescentPrev.sNew(c,:);
    runDescentOnTheta =  DescentPrev.runDescentOnTheta(c,:);
    sPrev =  DescentPrev.sCurrent(c,:);
    fdCostPrev = DescentPrev.fdCost{c};
    descentStep = DescentPrev.descentStep{c};
    iRepeat = DescentPrev.iRepeatNext;
    iDisc = DescentPrev.iDisc;
    

    
end



[fdCost] = FiniteDifferencePositiveGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M);
% % load('temp2.mat')
% vecnorm(fdCost,2,1)

iStop = find(abs(vecnorm(fdCost,2,1))<1E-2 | abs(vecnorm(fdCost,2,1))>1E100); %terminate continuation of descent
runDescentOnTheta(iStop) = 0;
% 


thetaNew = thetaCurrent;
igradientTheta = runDescentOnTheta & ~iDisc;
thetaNew(:,igradientTheta) =  thetaCurrent(:,igradientTheta)-descentStep(:,igradientTheta).*fdCost(:,igradientTheta); %calculate new descent cost
iDiscTheta = runDescentOnTheta & iDisc;
thetaNew(:,iDiscTheta) =  thetaCurrent(:,iDiscTheta)-descentStep(:,iDiscTheta).*sign(fdCost(:,iDiscTheta)); %calculate new descent cost





% if DescentPrev.Count > 0 & ~DescentPrev.newStart 
%     [thetaNew,descentStep,runDescentOnTheta,iDisc] = AdaptiveStepCheck(theta,S,thetaCurrent,thetaNew,fdCost,descentStep,runDescentOnTheta,sPrev,iDisc,M);
% end

% figure(100)
% thetaNorm = vecnorm(theta,2,1);
% plot(thetaNorm,S,'.')
% hold on
% 
% plot(vecnorm(thetaCurrent,2,1),sCurrent,'o')


if any(runDescentOnTheta)
    %Run Rise and Fall Method for new initial conditions
    [sNew,phiSetNew] = CostFunction(thetaNew(:,runDescentOnTheta),M);
    
    %store output data
    thetaOut = [theta thetaNew(:,runDescentOnTheta)];
    phiSetOut = [phiSet phiSetNew];
    sOut = [S  sNew];
    sNewOut = sPrev; sNewOut(runDescentOnTheta) = sNew;
    
    thetaOutNorm = vecnorm(thetaOut,2,1);
    [thetaOutNorm,iSortTheta] = sort(thetaOutNorm,'ascend');
    [thetaOut] = thetaOut(:,iSortTheta);
    phiSetOut = phiSetOut(iSortTheta);
    sOut = sOut(iSortTheta);
    TerminateFlag = false;
else
    thetaOut = theta;
    thetaOutNorm = vecnorm(thetaOut,2,1);
    phiSetOut = phiSet;
    sOut = S;
    sNewOut = sPrev;
    TerminateFlag = true;
end

% plot(vecnorm(thetaNew,2,1),sNewOut,'x','color','k')
% drawnow()
% hold off


%Allow cancellation of bad moves
iCancelLastMove = (sNewOut>=1.2*sPrev); %indices that increased the cost by more than 20% 
thetaNew(:,iCancelLastMove) = thetaCurrent(:,iCancelLastMove);
sNewOut(iCancelLastMove) = sPrev(iCancelLastMove);
% fdCost(:,iCancelLastMove) = fdCostPrev(:,iCancelLastMove);
iRepeatNext = false(size(sCurrent));
iRepeatNext(iCancelLastMove) = true;
descentStep(:,iCancelLastMove) = descentStep(:,iCancelLastMove)/2;
for i = 1:size(descentStep,2)
    if any(descentStep(:,i)<M.descent.minGamma)
        runDescentOnTheta(i) = 0;
    end
end

%Allow step change given not too bad moves
iStepChange = (sNewOut>sPrev); %indices that increased the cost 
iStepChange = iStepChange & ~iCancelLastMove;
descentStep(:,iStepChange) = descentStep(:,iStepChange)/2;



if DescentPrev.newStart
    DescentPrev.Count = 0;
else
    Descent = DescentPrev;
end
Descent.newStart = false;
Descent.Count = DescentPrev.Count+1;
c = Descent.Count;

Descent.runDescentOnTheta(c,:) = runDescentOnTheta;
Descent.theta{c} = thetaOut;
Descent.S{c} = sOut;
Descent.thetaCurrent{c}= thetaCurrent;
Descent.sCurrent(c,:) = sCurrent;
Descent.thetaNew{c} = thetaNew;
Descent.sNew (c,:) = sNewOut;
Descent.fdCost{c} = fdCost;
Descent.descentStep{c} = descentStep;
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


% 
function [fdCost] = FiniteDifferencePositiveGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M)
    d = M.dim-1;
    fdStep = M.descent.fdStep*ones(size(sCurrent));
    fdStep(iDisc) = descentStep(iDisc);
    iFindFDCost = find(~iRepeat & runDescentOnTheta); %check which indices need new gradient
    fdCost = fdCostPrev; %load prior gradient values
    if any(iFindFDCost)
        fdTheta = zeros(d,d);
        fdThetaAll = zeros(d,d*length(iFindFDCost));
        fdSInit = zeros(1,d*length(iFindFDCost));
        fdStepDivisor =  zeros(1,d*length(iFindFDCost));
        for ii = 1:length(iFindFDCost)
            for id = 1:d
                iFD = iFindFDCost(ii);
                fdStepVector = zeros(d,1);
                fdStepVector(id,1) = fdStep(iFD);
                fdTheta(:,id) = thetaCurrent(:,iFD)+fdStepVector;
            end
            iA = 1+(ii-1)*d;
            iB = iA+(d-1);
            fdThetaAll(:,iA:iB) = fdTheta;
            fdSInit(1,iA:iB) = repmat([sCurrent(iFD)],1,d);
            fdStepDivisor(1,iA:iB) = repmat([fdStep(iFD)],1,d);
        end 
        [fdS,~] = CostFunction(fdThetaAll,M); %cost near initial point 
        fdGradient = (fdS - fdSInit)./fdStepDivisor;
        fdCost(:,iFindFDCost) = reshape(fdGradient,d,length(iFindFDCost));
    end
end

function [thetaNew,descentStep,runDescentOnTheta,iDisc] = AdaptiveStepCheck(theta,S,thetaCurrent,thetaNew,fdCost,descentStep,runDescentOnTheta,sPrev,iDisc,M)

 
iCheck = SmartChecks3Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sPrev);
    iDisc = iDisc | iCheck;
    %adaptive step to account for smart checks
    if any(iCheck)

       while(any(iCheck))
           iNeedsAdapt = find(iCheck);
           for j = 1:length(iNeedsAdapt)
               jAdapt = iNeedsAdapt(j);
               if all(descentStep(:,jAdapt) > M.descent.minGamma) %limit to the adaptive step, if it shrinks too much end this descent
                 descentStep(:,jAdapt) = descentStep(:,jAdapt)./1.01;
                 thetaNew(:,jAdapt) =  thetaCurrent(:,jAdapt)-descentStep(:,jAdapt).*sign(fdCost(:,jAdapt));
               else
                   fprintf("Terminating due to small adaptive step\n")
                   thetaNew(:,jAdapt) = thetaCurrent(:,jAdapt);
                   runDescentOnTheta(jAdapt) = 0;
               end
           end
           iCheck = SmartChecks3Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sPrev);
       end
%        fprintf("Using an adaptive step, smallest step = %e\n", min(descentStep))
    end
    if any(ismember(thetaNew(:,runDescentOnTheta)',theta','rows'))
            error("New theta is the same as old theta\n")
    end
end


function [SNew,phiSetNew] = CostFunction(theta,M)
%     M.progressbar = false;
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

function [iCheck] = SmartChecks3Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sOld)

%         SPredict = griddatan(theta',S',thetaNew','nearest');
        
        iCheck2 = ismember(thetaNew',theta','rows')';

%         iCheck = SPredict'>=1.2*sOld;
        iCheck = (iCheck2) & runDescentOnTheta ;


end

