function [M] = Parameters(parameterNames,parameterValues)
%Configure PARAMETERS here 
note = "Two Forced Duffing Oscillators";
paramNote = "Two Oscillator,Theta GridSearch";
a1 = 1; a3 = .3; nu = .1; F = .4; w = 1.4; kc = .1; %rhs parameters (note basin interpolant mat file must be changed if rhs parameters are changed)
dim = 4; %deterministic system dimension
rIC = 10^-15; %radius of momenta initial conditions
pp = 0; %poincare phase
nIC = 10; %number of initial conditions 
rhsString = 'TwoDuffing';
 T = 2*pi/w;  dT = T/2; dt = T/32; tf = 500*T;
solver = @ode45;
psiEps = .05; %phase threshold
tFall = 10*T; %Max amount of time for system to fall to attractor radius
plotFall = true;
fastPostProcessing = false;
rA1 = .3; %radius of initial sphere around initial attractor
rA = .1; %accepted radius around an attractor
rS =.001; %accepted radius around a saddle
tstep = .1; %time that must pass prior to checking for sphere condition
iA = 3; %initial attractor fixed point identifier
onceAPeriod = true;
terminateType = 'DuffingBoundary'; 
nWorkers = Inf;
continueRun  = false;
clusterRun = CheckIfCluster();
xcoordinates = true;
uniformInX = true;
nRVs = 2000; %number random variables per dimension for random IC initialization
saveMemory = 1;
methodTest = false;
searchAlgorithm = "Gradient Descent Discontinuous 2";



%One oscillator modications
% iA = 1;
% rhsString = 'Duffing';
% dim = 2; 
% note = "One Forced Duffing Oscillator";
% paramNote = "One Oscillator";
% nIC = 10;

% Test modifications
methodTest = true;
saveMemory = true;
dim = 2;
note = "methodTest";
paramNote = "methodTest";
nIC = 10;
xcoordinates = false;



%MinSearch Parameters
nLM = 10; %maximum number of local minimum to explore
maxIter = 40;

%Descent parameters
if contains(searchAlgorithm,"Gradient") || contains(searchAlgorithm,"Fletcher")
    descent.Gamma = .1; 
else
    descent.Gamma = 1; 
end
descent.fdStep = 1E-4; %finite difference step
descent.minGamma = 1E-10;
descent.discGamma = .5;
descent.DiscThresh = 2.0;
progressbar = true;

descent.stochasticNonGradientSearch = true;



%Calculated parameters
% D = 2*length(qo);
D = 2*dim;
tspan = [0:dt:tf];


%Store Parameters in a structure
M.note = note;
M.rIC = rIC; 
M.nIC = nIC;
M.rhsString = rhsString;
M.tf = tf; M.dt = dt;
M.solver = solver;
M.D = D;
M.tspan = tspan;
% M.theta = theta; M.dtheta = dtheta;
M.dim = dim;
M.terminateType = terminateType;
M.nWorkers = nWorkers;
M.plotFall = plotFall;
M.fastPostProcessing = fastPostProcessing;
M.pp = pp;
M.dT = dT;
M.progressbar = progressbar;
M.continueRun = continueRun;
M.clusterRun = clusterRun;
M.xcoordinates = xcoordinates;
M.uniformInX = uniformInX;
M.transMatrix = [];
M.nRVs = nRVs;
M.saveMemory = saveMemory;
M.methodTest = methodTest;
M.searchAlgorithm = searchAlgorithm;


M.paramNote = paramNote;
M.Mrhs.a1 = a1; M.Mrhs.a3 = a3; M.Mrhs.nu = nu; M.Mrhs.kc = kc;
M.Mrhs.F = F; M.Mrhs.w = w; M.Mrhs.psiEps = psiEps;
 M.Mrhs.rA = rA; M.Mrhs.rS = rS; M.Mrhs.tstep = tstep; %M.Mrhs.qo = qo;
M.Mrhs.rA1 = rA1; M.Mrhs.tFall = tFall; M.Mrhs.iA = iA; M.Mrhs.T = T;
M.Mrhs.onceAPeriod = onceAPeriod; M.Mrhs.dim = dim;

M.MS.nLM = nLM; M.MS.maxIter = maxIter;

M.descent = descent;


%modify parameter values via function input during runtime
if nargin == 2
    for i = 1:length(parameterNames)
        if contains(parameterNames{i},'MRHS_')
            PN = strrep(parameterNames{i},'MRHS_','');
            M.Mrhs = setfield(M.Mrhs,PN,parameterValues{i});
        else
            M = setfield(M,parameterNames{i},parameterValues{i});
        end
    end
end


%Set Function Handles
M.RHS = str2func([M.rhsString,'RHS']);
M.HamiltonianRHS = str2func([M.rhsString,'HamiltonianRHS']);
M.Lagrangian = str2func([M.rhsString,'Lagrangian'])
% M.TerminateEvent = str2func(['TerminateAt',M.terminateType,'Event']);

M.Mrhs.RHS = M.RHS; M.Mrhs.solver= solver;
end

