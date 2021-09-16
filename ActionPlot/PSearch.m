function [phiSetOut,msLog,TerminateFlag] = PSearch(phiSet,msLog,nLM,M)
%MINSEARCH1D                 

theta = msLog{end}{1};
S = msLog{end}{2};
DescentPrev = msLog{end}{3};


[thetaCurrent,sCurrent,~] = IdentifyLocalMinima(theta,S,nLM);

runDescentOnTheta = true(size(thetaCurrent));
descentStep = M.descent.Gamma*ones(size(thetaCurrent));
sPrev = sCurrent;
fdCostPrev = NaN(size(thetaCurrent));
iRepeat = false(size(thetaCurrent));
options = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
lb = [0];
ub = [2*pi];
A = [];
b = [];
Aeq = [];
beq = [];
fun = @(theta) CostFunction(theta,M);
nonlcon = @mycon;
for i = 3;
    i
    x0 = thetaCurrent(i);
    [x,fval,exitflag,output] = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

phiSetOut = phiSet;
TerminateFlag = true;
end 

function [SNew] = CostFunction(theta,M)
    M.progressbar = false;
     phiSetNew = PostProcessTrajectories2(IntegrateRHS(GenerateInitialConditionsFloquet(theta,M),M),M);
     SNew = IntegrateLagrangian(phiSetNew,M); 
end
function [c,ceq] = mycon(x)
 c= [];
 ceq = [];
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