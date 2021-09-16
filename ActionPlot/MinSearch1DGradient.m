function [phiSetOut,msLog] = MinSearch1D(phiSet,msLog,nLM,M)
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

%Neighbor Differentials
ds1 = [0 ;diff(sAug)]; 
ds2 = [diff(sAug);0];


iLM = find(ds1<0 & ds2>0); %index of local minima
sLM = S(iLM); %energy of local minima

%Limit to max number of local minima to use
if length(sLM) > nLM
    [~,iSort] = sort(sLM,'ascend');
    iLM = iLM(iSort);
    iLM = iLM(1:nLM);
    sLM = S(iLM);
end


subgamma = .01;
gamma = .1;

thetai = thetaAug(iLM);

Sright = CostFunction(thetai+subgamma,M);
FD  = (Sright-sLM)/subgamma

iStop = find(abs(FD)<5E-2 | abs(FD)>1E1)
thetaStop = thetai(iStop);

thetai(iStop) = []; FD_Copy = FD; FD(iStop) = [];

thetaNew =  thetai-gamma*FD';

iRepeat = ismember(thetaNew,thetaAug);

%break up line into 8ths near discontinuity
gamma_line8th = gamma./8;
if any(iRepeat)
   fprintf("Drawing line near possible discontinuity\n")
   thetaStart = thetai(iRepeat);
   fdStart(1,:) = FD(iRepeat);
   thetaNew(iRepeat) = [];
   for i = 1:length(thetaStart)
      for j = 1:7
         thetaNew(end+1) =  thetaStart(i)-j*gamma_line8th*fdStart(i);
      end
   end
end



% %generate new initial conditions
% iLeft = iLM - 1; iRight = iLM +1;
% 
% thetaNew =  mean([thetaAug(iLM) thetaAug(iLM); thetaAug(iLeft) thetaAug(iRight)]);

%Run Rise and Fall Method for new initial conditions
[phiSetNew,SNew] = GetNewPhiSets(thetaNew,M);

%store output data
thetaOut = [theta thetaNew];
phiSetOut = [phiSet phiSetNew];
Sout = [S ; SNew];

[thetaOut,iSortTheta] = sort(thetaOut,'ascend');
phiSetOut = phiSetOut(iSortTheta);
Sout = Sout(iSortTheta);

Descent.iLocalMinima = find(ismember(thetaOut,thetaAug(iLM)));
Descent.thetaMinima = thetaOut(Descent.iLocalMinima);
Descent.FD = FD_Copy';
Descent.thetaStop = [DescentPrev.thetaStop thetaStop];;


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

