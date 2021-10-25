function [phiSetOut,msLog,TerminateFlag] = MinSearchGradientDescent(phiSet,msLog,nLM,M)
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

%     [~,sCurrent,iLM] = IdentifyLocalMinima(theta,S,nLM);
    [~,sCurrent,iLM] = IdentifySmallValues(theta,S,nLM);
    thetaCurrent = theta(:,iLM);
    runDescentOnTheta = true(size(sCurrent));
    descentStep = M.descent.Gamma*ones(size(thetaCurrent));
    sPrev = sCurrent;
    fdCostPrev = NaN(size(thetaCurrent));
    iRepeat = false(size(sCurrent));
    iDisc = false(size(sCurrent));
    thetaDescentIndex = zeros(size(S));

    
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
%     thetaDescentIndex = DescentPrev.thetaDescentIndex;
    

    
end

%Descent Computation
if M.xcoordinates
    [fdCost] = FiniteDifferencePositiveGradientOnXCoordinates(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M);
else
    [fdCost] = FiniteDifferencePositiveGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M);
end

iStop = find(abs(vecnorm(fdCost,2,1))<1E-2);
iDiscontinuity = find(abs(vecnorm(fdCost,2,1))>1E2); %terminate continuation of descent
runDescentOnTheta(iStop) = false;
iDisc(iDiscontinuity) = true;
% 



thetaNew = thetaCurrent;
igradientTheta = runDescentOnTheta & ~iDisc;
thetaNewGradient =  thetaCurrent(:,igradientTheta)-descentStep(:,igradientTheta).*fdCost(:,igradientTheta); %calculate new descent cost
iDiscTheta = runDescentOnTheta & iDisc;
% % thetaNewGridSearch =  thetaCurrent(:,iDiscTheta)-descentStep(:,iDiscTheta).*sign(fdCost(:,iDiscTheta)); %calculate new descent cost


[thetaNewGridSearch,iGridSearch] = GridSearchNearDiscontinuity(theta,S,thetaCurrent(:,iDiscTheta),find(iDiscTheta),sCurrent(iDiscTheta),M);

% if DescentPrev.Count > 0 & ~DescentPrev.newStart 
%     [thetaNew,descentStep,runDescentOnTheta,iDisc] = AdaptiveStepCheck(theta,S,thetaCurrent,thetaNew,fdCost,descentStep,runDescentOnTheta,sPrev,iDisc,M);
% end

if M.xcoordinates
    thetaNew = thetaNew./vecnorm(thetaNew,2,1);
else
    thetaNew = mod(thetaNew,2*pi);
    thetaNew(2:3,:) = mod(thetaNew(2:3,:),pi);
end

% if ~CheckIfCluster()
%     figure(100)
%     thetaNorm34 = vecnorm(theta(3:4,:),2,1);
%     scatter3(theta(1,:),theta(2,:),sign(theta(4,:)).*thetaNorm34,40,S,'filled')
%     cb = colorbar;
%     hold on
%     plot(vecnorm(thetaCurrent,2,1),sCurrent,'o')
% end




if any(igradientTheta) | any(iDiscTheta)
    
    sNewOut = sPrev;
    sOut = S;
    phiSetOut = phiSet;
    thetaOut = theta;
    if any(igradientTheta)
        %Run Rise and Fall Method for new initial conditions
        [sNew,phiSetNew] = CostFunction(thetaNewGradient,M);

        %store output data
        thetaNew(:,igradientTheta) = thetaNewGradient;
        sNewOut(igradientTheta) = sNew;
        [thetaOut,phiSetOut,sOut] = MergeNewData(thetaOut,phiSetOut,sOut,thetaNewGradient,phiSetNew,sNew);
        

         
    end
    
    if any(iDiscTheta)
        [sNew,phiSetNew] = CostFunction(thetaNewGridSearch,M);
        
        iGridSearchUnique = unique(iGridSearch);
        for i = 1:length(iGridSearchUnique)
            ii = find(iGridSearchUnique(i) == iGridSearch);
            [~,imin] = min(sNew(ii));
            iExtendedDiscTheta(i) = ii(imin);
        end
         thetaNew(:,iDiscTheta) = thetaNewGridSearch(:,iExtendedDiscTheta);
        sNewOut(iDiscTheta) = sNew(iExtendedDiscTheta);
        [thetaOut,phiSetOut,sOut] = MergeNewData(thetaOut,phiSetOut,sOut,thetaNewGridSearch,phiSetNew,sNew);
    end
    
    

    TerminateFlag = false;
else
    thetaOut = theta;
    phiSetOut = phiSet;
    sOut = S;
    sNewOut = sPrev;
    TerminateFlag = true;
end

% if ~CheckIfCluster()
%     scatter3(thetaNew(1,:),thetaNew(2,:),sign(thetaNew(4,:)).*vecnorm(thetaNew(3:4,:),2,1),'xr')
%     drawnow()
%     hold off
% end


%Allow cancellation of bad moves
iCancelLastMove = (sNewOut>=M.descent.DiscThresh*sPrev); %indices that increased the cost by more than 20% 
thetaNew(:,iCancelLastMove) = thetaCurrent(:,iCancelLastMove);
sNewOut(iCancelLastMove) = sPrev(iCancelLastMove);
% fdCost(:,iCancelLastMove) = fdCostPrev(:,iCancelLastMove);
iRepeatNext = false(size(sCurrent));
iRepeatNext(iCancelLastMove) = true;
descentStep(:,iCancelLastMove) = descentStep(:,iCancelLastMove)/2;
iDisc(iCancelLastMove) = true;
for i = 1:size(descentStep,2)
    if any(descentStep(:,i)<M.descent.minGamma)
        runDescentOnTheta(i) = false;
%     elseif any(descentStep(:,i)<M.descent.discGamma)
%         iDisc(i) = true
    end
end

%Allow step change given not too bad moves
iStepChange = (sNewOut>sPrev); %indices that increased the cost 
iStepChange = iStepChange & ~iCancelLastMove;
descentStep(:,iStepChange) = descentStep(:,iStepChange)/2;

descentStep

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
% Descent.thetaDescentIndex = thetaDescentIndex;
    

msLog{end+1} = {thetaOut,sOut,Descent};

% vecnorm(thetaNew,2,1)
sCurrent
sNewOut
end

function [thetaOut,phiSetOut,sOut] = MergeNewData(theta,phiSet,S,thetaNew,phiSetNew,sNew)
    thetaOut = [theta thetaNew];
    phiSetOut = [phiSet phiSetNew];
    sOut = [S  sNew];
    
%     thetaOutNorm = vecnorm(thetaOut,2,1);
%     [thetaOutNorm,iSortTheta] = sort(thetaOutNorm,'ascend');
%     [thetaOut] = thetaOut(:,iSortTheta);
%     phiSetOut = phiSetOut(iSortTheta);
%     sOut = sOut(iSortTheta);

end

function [thetaOut,sOut,iOut] = IdentifySmallValues(theta,S,nLM)

% %Augmented optimization space
%     thetaAug = [theta theta(1)];
%     sAug = [S S(1)]; %adding the first term to the end to complete the circle
% 
%  %Find local minima
%     %Neighbor Differentials
%     ds1 = [0 diff(sAug)]; 
%     ds2 = [diff(sAug) 0];
% 
% 
%     iLM = find(ds1<0 & ds2>0); %index of local minima
%     sLM = S(iLM); %energy of local minima
% 
%     %Limit to max number of local minima to use
%     if length(sLM) > nLM
%         [~,iSort] = sort(sLM,'ascend');
%         iLM = iLM(iSort);
%         iLM = iLM(1:nLM);
%         sLM = S(iLM);
%     end
%     thetaLM = thetaAug(iLM);

    
    [~,is] = sort(S,'ascend');
    if length(is) > nLM
        is = is(1:nLM);
    end
    thetaOut = theta(:,is);
    sOut = S(is);
    iOut = is;
end

function [thetaOut,sOut,iOut] = IdentifyLocalMinima(theta,S,nLM)


    thetaAug = [theta];
    sAug = [S]; 
   
    [s,is] = sort(S','ascend');
    x = theta(:,is)';
    
    ix = 1; n = 0;
    eps = .01; 
    eps_perm = [eye(4);-eye(4)]*eps;
    while(n<=nLM & ix<size(theta,2))
        xi = x(ix,:);
        si = s(ix);
        for i = 1:size(eps_perm,1)
            xq(i,:) = xi+eps_perm(i,:);
        end
        sq = griddatan(x,s,xq)
        
        
    end
    
   dummy = 1; 
% 
%  %Find local minima
%     %Neighbor Differentials
%     ds1 = [0 diff(sAug)]; 
%     ds2 = [diff(sAug) 0];
% 
% 
%     iLM = find(ds1<0 & ds2>0); %index of local minima
%     sLM = S(iLM); %energy of local minima
% 
%     %Limit to max number of local minima to use
%     if length(sLM) > nLM
%         [~,iSort] = sort(sLM,'ascend');
%         iLM = iLM(iSort);
%         iLM = iLM(1:nLM);
%         sLM = S(iLM);
%     end
%     thetaLM = thetaAug(iLM);

    
end
%
function [fdCost] = FiniteDifferencePositiveGradient(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M)
    d = M.dim-1;
    fdStep = M.descent.fdStep*ones(size(sCurrent));
    fdStep(iDisc) = descentStep(iDisc);
    iFindFDCost = find(~iRepeat & runDescentOnTheta &~iDisc); %check which indices need new gradient
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


function [fdCost] = FiniteDifferencePositiveGradientOnXCoordinates(thetaCurrent,sCurrent,runDescentOnTheta,fdCostPrev,iRepeat,iDisc,descentStep,M)
    % x1^2+x2^2+x3^2+x4^2 = 1;
    % df1 = 2x1; df2 = 2x2; df3 = 2x3; df4 = 2x4;
    % p = point, y = coordinate
    % A = df/norm(df); %radial vector
    % B = null(A)'; orthogonal vectors to radial
    % NewX = X + fdStep*gradS
    % NewTheta = ConverToTheta(X)
    
    

    d = M.dim;
    dT = d - 1; %d tangent space
    
    
    fdStep = M.descent.fdStep*ones(size(sCurrent));
    fdStep(iDisc) = descentStep(iDisc);
    iFindFDCost = find(~iRepeat & runDescentOnTheta & ~iDisc); %check which indices need new gradient
    fdCost = fdCostPrev; %load prior gradient values
    if any(iFindFDCost)
        fdTheta = zeros(d,d-1);
        fdThetaAll = zeros(d,dT*length(iFindFDCost));
        fdSInit = zeros(1,dT*length(iFindFDCost));
        fdStepDivisor =  zeros(1,dT*length(iFindFDCost));
        for ii = 1:length(iFindFDCost)
            iFD = iFindFDCost(ii);
            P = thetaCurrent(:,iFD);
            eRadial = (2*P)./norm(P);
            eOrthogonal = null(eRadial')';
            for id = 1:dT %iterate tangent space dimension
                fdStepVector(:,id) =  fdStep(iFD)*eOrthogonal(id,:)';
                fdTheta(:,id) = thetaCurrent(:,iFD)+fdStepVector(:,id);
            end
            iA = 1+(ii-1)*dT;
            iB = iA+(dT-1);
            fdThetaAll(:,iA:iB) = fdTheta;
            fdSInit(1,iA:iB) = repmat([sCurrent(iFD)],1,dT);
            fdOrthogonalVectors(:,iA:iB) = eOrthogonal';
            fdStepDivisor(1,iA:iB) = fdStep(iFD);
        end 
        fdThetaAllNorm = vecnorm(fdThetaAll,2,1);
        fdThetaAllNorm = repmat(fdThetaAllNorm,d,1);
        fdThetaAll = fdThetaAll./fdThetaAllNorm;
        [fdS,~] = CostFunction(fdThetaAll,M); %cost near initial point 
        fdGradientExpanded = repmat((fdS - fdSInit)./fdStepDivisor,d,1).*fdOrthogonalVectors;
        for ii = 1:length(iFindFDCost);
            iA = 1+(ii-1)*dT; iA2 = (ii)*dT;
            fdGradient(:,ii) = sum(fdGradientExpanded(:,iA:iA2),2);
        end
        fdCost(:,iFindFDCost) = fdGradient;
    end
end

function [thetaNew,iGridSearch] = GridSearchNearDiscontinuity(theta,s,thetaNearDisc,iDiscThetaFind,sNearDisc,M)
    d = size(thetaNearDisc,1);
    nNewPoints = d*2-1;
    thetaNew = zeros(d,nNewPoints*size(thetaNearDisc,2));
    iGridSearch = zeros(1,nNewPoints*size(thetaNearDisc,2));
    for i = 1:size(thetaNearDisc,2)
        iNewStart = (i-1)*(nNewPoints)+1;
        distance = vecnorm(theta-thetaNearDisc(:,i),2,1);
        [distanceSorted,iSortDistance] = sort(distance,'ascend');
        sSorted = s(iSortDistance);
        thetaSorted = theta(:,iSortDistance);
        iClosestDiscPoint = min(find(sSorted>M.descent.DiscThresh*sNearDisc(i)));
        if isempty(iClosestDiscPoint)
            continue
        end
        thetaAtNearDisc = thetaSorted(:,iClosestDiscPoint);
        dirMagnitude = distanceSorted(iClosestDiscPoint);
        dirVector = (thetaAtNearDisc-thetaNearDisc(:,i))./dirMagnitude;
        orthogonalVectors = null(dirVector');
        L = size(orthogonalVectors,2);
        
        thetaNew(:,iNewStart) = thetaNearDisc(:,i)+.5*dirMagnitude*dirVector;
        thetaNew(:,iNewStart+1:iNewStart+L) = thetaNew(:,iNewStart)+.1*dirMagnitude*orthogonalVectors;
        thetaNew(:,iNewStart+1+L:iNewStart+2*L) = thetaNew(:,iNewStart)-.1*dirMagnitude*orthogonalVectors;
        iGridSearch(iNewStart:iNewStart+2*L) = iDiscThetaFind(i);
    end
        iGridSearch(iGridSearch==0) = [];
        thetaNew(:,all(thetaNew==0,1)) = [];
       
        
    
end

function [thetaNew,descentStep,runDescentOnTheta,iDisc] = AdaptiveStepCheck(theta,S,thetaCurrent,thetaNew,fdCost,descentStep,runDescentOnTheta,sPrev,iDisc,M)

 
iCheck = SmartChecks3Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sPrev,thetaCurrent);
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
           iCheck = SmartChecks3Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sPrev,thetaCurrent);
       end
%        fprintf("Using an adaptive step, smallest step = %e\n", min(descentStep))
    end
    if any(ismember(thetaNew(:,runDescentOnTheta)',theta','rows'))
            error("New theta is the same as old theta\n")
    end
end


function [iCheck] = SmartChecks1Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sOld)

        for i = 1:length(thetaNew)
            SPredict(1,i) = interp1(theta,S,thetaNew(i),'nearest');
        end
        iCheck2 = ismember(thetaNew,theta);

        iCheck = SPredict>=M.descent.DiscThresh*sOld;
        iCheck = (iCheck | iCheck2) & runDescentOnTheta ;


end

function [iBadStep] = SmartChecks3Ddiscontinuity(theta,S,thetaNew,runDescentOnTheta,sOld,thetaCurrent)

        SPredict = griddatan(theta',S',thetaNew','nearest');
        iCheck = SPredict'>=M.descent.DiscThresh*sOld;
        if any(iCheck)
            iThetaNearest = knnsearch(theta',thetaNew(:,iCheck)');
            dCurrentToNewThetaNearestNeighbor = vecnorm(theta(:,iThetaNearest)-thetaCurrent(:,iCheck),2,1);
            dCurrentToNewTheta = vecnorm(thetaNew(:,iCheck)-thetaCurrent(:,iCheck),2,1);
            iCheck2 = dCurrentToNewThetaNearestNeighbor > dCurrentToNewTheta;
            iCheck(iCheck) = iCheck2;
        end
        
        iCheck3 = ismember(thetaNew',theta','rows')';

%         
        iBadStep = (iCheck | iCheck3) & runDescentOnTheta ;


end

