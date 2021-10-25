function [phiSetOut,msLog,TerminateFlag] = MinSearchGridSearch(phiSet,msLog,nLM,M)
%MINSEARCH1D                 
% Identify all local minima (all nearest neighbors have larger energy)
% Select up to nLM local minima with smallest energy
% Generate initial conditions between the neighbors to the local minima and the local minima.
% Run rise and fall method for new initial conditions


%extract data from min search log
theta = msLog{end}{1};
s = msLog{end}{2};
gsPrev = msLog{end}{3};




    if gsPrev.newStart  %first loop

        [~,sCurrent,iLM] = IdentifySmallValues(theta,s,nLM);
        thetaCurrent = theta(:,iLM);
        runGSOnTheta = true(size(sCurrent));
        stepSize = M.descent.Gamma*ones(size(thetaCurrent));
    else
        c = gsPrev.Count;
        thetaCurrent = gsPrev.thetaNew{c};
        sCurrent =  gsPrev.sNew(c,:);
        runGSOnTheta =  gsPrev.runGSOnTheta(c,:);
        stepSize = gsPrev.stepSize{c};

    end
%     PL_ThetaSpace2(theta,s)
%     PL_ThetaSpace3(thetaCurrent)
%     view(0,0)
%     drawnow()
%     hold off
%     dummy = 1;
    if M.xcoordinates
        [thetaNewGridSearch,iGridSearch] = GridSearchOnX(thetaCurrent(:,runGSOnTheta),find(runGSOnTheta),stepSize,M);
    else
        [thetaNewGridSearch,iGridSearch] = GridSearch(thetaCurrent(:,runGSOnTheta),find(runGSOnTheta),stepSize,M);
        [thetaNewGridSearch] = ModTheta(thetaNewGridSearch)
    end



    
    [sNewGridSearch,phiSetNewGridSearch] = CostFunction(thetaNewGridSearch,M);


    iGridSearchUnique = unique(iGridSearch);
    iNewValue = zeros(1,length(iGridSearchUnique));
    for i = 1:length(iGridSearchUnique)
        ii = find(iGridSearchUnique(i) == iGridSearch);
        [~,imin] = min(sNewGridSearch(ii));
        iNewValue(i) = ii(imin);
    end
    
    thetaNew = thetaCurrent;
    sNew = sCurrent;
    thetaNew(:,runGSOnTheta) = thetaNewGridSearch(:,iNewValue);
    sNew(runGSOnTheta) = sNewGridSearch(iNewValue);
    [thetaOut,phiSetOut,sOut] = MergeNewData(theta,phiSet,s,thetaNewGridSearch,phiSetNewGridSearch,sNewGridSearch);
   


    

    iCancelMove = (sNew>=sCurrent);
    thetaNew(:,iCancelMove) = thetaCurrent(:,iCancelMove);
    sNew(iCancelMove) = sCurrent(iCancelMove);
    stepSize(iCancelMove) = stepSize(iCancelMove)/2;
    
%     [maxS] = max(sNew(runGSOnTheta)); %each iteration kill worst solution branch
%     iEndBranch = find(maxS == sNew);
%     if length(iEndBranch) == 1
%         runGSOnTheta(iEndBranch) = false;
%     end
    
    iStop = find(stepSize<= 1E-5);
    runGSOnTheta(iStop) = false;
    
    TerminateFlag = ~any(runGSOnTheta);
    

    if gsPrev.newStart
        gsPrev.Count = 0;
    else
        gs = gsPrev;
    end
    gs.newStart = false;
    gs.Count = gsPrev.Count+1;
    c = gs.Count;

    gs.runGSOnTheta(c,:) = runGSOnTheta;
    gs.theta{c} = thetaOut;
    gs.S{c} = sOut;
    gs.thetaCurrent{c}= thetaCurrent;
    gs.sCurrent(c,:) = sCurrent;
    gs.thetaNew{c} = thetaNew;
    gs.sNew (c,:) = sNew;
    gs.stepSize{c} = stepSize;
    msLog{end+1} = {thetaOut,sOut,gs};

end

function [thetaOut,phiSetOut,sOut] = MergeNewData(theta,phiSet,S,thetaNew,phiSetNew,sNew)
    thetaOut = [theta thetaNew];
    phiSetOut = [phiSet phiSetNew];
    sOut = [S  sNew];

end

function [thetaOut,sOut,iOut] = IdentifySmallValues(theta,S,nLM)

    [~,is] = sort(S,'ascend');
    if length(is) > nLM
        is = is(1:nLM);
    end
    thetaOut = theta(:,is);
    sOut = S(is);
    iOut = is;
end

%

function [thetaNew,iGridSearch] = GridSearch(thetaCurrent,iCurrent,stepSize,M)
   
    d = size(thetaCurrent,1);
    n = size(thetaCurrent,2)
    nNewPoints = d*2;
    thetaNew = zeros(d,nNewPoints*n);
    iGridSearch = zeros(1,nNewPoints*n);
   
    for i = 1:n
        j = (i-1)*(nNewPoints)+1;
        k = (i)*(nNewPoints);
        if M.descent.stochasticGridSearch
            v1 = rand(1,d);
        else
            v1 = zeros(1,d); v1(1) = 1;
        end
        v = [v1'/norm(v1) null(v1)];
        step = stepSize(i)*[v -v];
        thetaNew(:,j:k) = thetaCurrent(:,i)+step;
        iGridSearch(j:k) = iCurrent(i);
    end
end

function [xNew,iGridSearch] = GridSearchOnX(xCurrent,iCurrent,stepSize,M)
   
    d = size(xCurrent,1);
    n = size(xCurrent,2)
    nNewPoints = (d-1)*2;
    xNew = zeros(d,nNewPoints*n);
    iGridSearch = zeros(1,nNewPoints*n);
   
    for i = 1:n
        j = (i-1)*(nNewPoints)+1;
        k = (i)*(nNewPoints);
        p = xCurrent(:,i)';
        vr = (2*p)./norm(p); %radial vector
        vt = null(vr); %tangent vectors
        c = repmat(rand(1,size(vt,2)),size(vt,1),1); 
        v1 = sum(c.*vt,2); v1 = v1'./norm(v1); %random tangent vector
        if M.descent.stochasticGridSearch
            v = [v1' null([vr;v1])]; %orthogonal vectors including random tangent vector
        else
            v = [vr' vt]; %orthogonal vectors not stochastic
        end
        step = stepSize(i)*[v -v];
        xNew(:,j:k) = xCurrent(:,i)+step;
        xNew(:,j:k) = xNew(:,j:k)./vecnorm(xNew(:,j:k),2,1);
        iGridSearch(j:k) = iCurrent(i);
    end
end


 

