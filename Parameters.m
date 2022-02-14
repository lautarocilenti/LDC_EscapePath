function [M] = Parameters(parameterNames,parameterValues) 
% Configure PARAMETERS here. 
% 02/11/2022 

% *************************************************************************
% Parameters
% *************************************************************************
% Most parameters can be kept as they are. To change system dimension see examples at the end of this file.
% Do note, attractors have been computed offline and saved as Matfiles that are loaded depending on paramters.
% If you want to run a  different system not provided you will have to compute your own
% attractors or modify the code to compute them during runtime.

note = "";
paramNote = "";

%deterministic system
a1 = 1; a3 = .3; nu = .1; F = .4; w = 1.4; kc = .01; %rhs parameters
dim = 2; %deterministic system dimension
rhsString = 'Duffing'; %Options {'Duffing','TwoDuffing','NDuffing'}

%Radius
rIC = 10^-10; %gamma - radius of Lagrangian Manifold Approximation
rA1 = .1; %radius of sphere around initial attractor (used to save comp. time by avoiding checks for basin of attraction very close to initial attractor)
rA = .1; %accepted radius around a new attractor
rS = 1E-5; %accepted radius around a saddle

%Attractors
iAString = "High Attractor"; %initial attractor fixed point identifier must match system dimension
fAString = "Any"; 
pp = 0; %theta_0 initial phase angle from 0 to 2*pi, only used when includePhase parameter is false 


%Optimization Parameters
nIC = 30; %number of initial conditions 
solver = @ode45; %Numerical integrator
searchAlgorithm = "Stochastic Grid"; %optimization solver Options ("NelderMead Simplex","Stochastic Grid","Fixed Grid","Gradienct Descent","Fletcher Reevers")
xcoordinates = false; %when true uses x coordinates for sphere in R^d. when false it uses angle coordinates for sphere in R^d. 'NelderMead Simplex' solver will not work when this is true 
uniformInX = true; %initial points are selected uniformly in the sphere. when false, initial points are selected uniformly in the angle space of the sphere
nRVs = 500; %number random variables in R^d for random IC initialization. Recommend range 100 - 15000 and always >nIC. Actual nIC initial conditions are downsampled uniformly from this initial set.  If the dimension is very large this will create an nRV X d matrix so be careful not to make too large.  
minICStopRemoval = 5; %optimization will run on nLM candidates for local minima in parallel. To reduce comp. time, it will throw out one candidate with highest cost each iteration but it will always keep at least minICStopRemoval 
nLM = 4; %maximum number of local minimum to explore (recommended 4 for testing, 30 to 50 for actual running) candidates are selected based on the initial guesses with lowest costs that are not too close to each other
maxIter = 3; %Stop Condition 1: maximum number of iterations of the optimizer. 
fdStep = 1E-2; %finite difference step size for approximating gradient. Too small or too large leads to errors in gradient estimate
minGamma = 1E-5; %Stop Condition 2: optimizer step size decreases below min gamma
discGamma = .5; %gradient descent switches to stochastic grid vector method after reaching discontinuity and resets the step size to this.
DiscThresh = 2.0; %cost multiplier to identify small step led to a discontinuous jump
includePhase = true; %increase the dimension of the optimization problem to include initial phase. true by default

%Plot Paramteters
plotFall = true; %includes the descent into a new attractor after laeaving initial basin

%Other parameters
continueRun  = false; %used to continue running optimization with known cost function values from previous run. Use with caution, will not recognize parameter changes to the system and two cost functions could be mixed in together.
clusterRun = CheckIfCluster(); %used to restricting plotting in enviroments that do not have GUI. Adapt the function to your method of detecting those enviroments
saveMemory = true; % when true, only stores initial points and not the entire trajectories
methodTest = false; %when true the cost function is replaced with a simple 2D cost function to test the solvers in 2D
costType = "action"; %Allows for modifying the cost function to evaluate other costs. Currently there are two options "distancefromSaddle","action". 'action' is always better bt the framework remains for testing other cost types
progressbar = true; %plots a progress bar in the text if true

%Calculated Parameters (paramters that are partially based on other parameters)
T = 2*pi/w; %Period
dT = T/2; %initial time between basin checks 
dt = T/32; %used to create time span but doesnt anything important since ode45 chooses the timestep of integration
tf = 500*T; %maximum integration (integration will terminate early if another basin is reached)
tFall = 10*T; %Basin check integration time interval
if contains(searchAlgorithm,"Gradient") || contains(searchAlgorithm,"Fletcher")
    Gamma = .1; %initial step size for optimizer
else
    Gamma = .25; %initial step size for optimizer
end
if contains(searchAlgorithm,"Gradient") || contains(searchAlgorithm,"Fletcher") || contains(searchAlgorithm,"Stochastic")
    stochasticNonGradientSearch = true;
else
     stochasticNonGradientSearch = false;
end
D = 2*dim; %hamiltonian system dimension
tspan = [0:dt:tf]; % time span for integration


% deprecated parameters  
psiEps = .005; %phase threshold
fastPostProcessing = false;
tstep = .1; %time that must pass prior to checking for sphere condition
onceAPeriod = true;
terminateType = 'DuffingBoundary'; 
nWorkers = Inf;
io = 0;










% *************************************************************************
% This section packages parameters into appropriate structures
% *************************************************************************

descent.Gamma  = Gamma;
descent.fdStep = fdStep; %finite difference step
descent.minGamma = minGamma;
descent.discGamma = discGamma;
descent.DiscThresh = DiscThresh;
descent.stochasticNonGradientSearch = stochasticNonGradientSearch;


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
M.costType = costType;
M.iAString = iAString;
M.fAString = fAString;
M.minICStopRemoval = minICStopRemoval;
M.includePhase = includePhase;


M.paramNote = paramNote;
M.Mrhs.a1 = a1; M.Mrhs.a3 = a3; M.Mrhs.nu = nu; M.Mrhs.kc = kc;
M.Mrhs.F = F; M.Mrhs.w = w; M.Mrhs.psiEps = psiEps;
M.Mrhs.rA = rA; M.Mrhs.rS = rS; M.Mrhs.tstep = tstep; %M.Mrhs.qo = qo;
M.Mrhs.rA1 = rA1; M.Mrhs.tFall = tFall; M.Mrhs.T = T;
M.Mrhs.onceAPeriod = onceAPeriod; M.Mrhs.dim = dim; M.io = io;

M.MS.nLM = nLM; M.MS.maxIter = maxIter;

M.descent = descent;

% *************************************************************************
% This section allows inputs to this function to modify parameter values during runtime. See the Superscripts folder for how this is used.
% *************************************************************************

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



% *************************************************************************
% The following are examples of modifications for different systems
% *************************************************************************


% One oscillator modications
% rhsString = 'Duffing';
% dim = 2; 
% note = "One Forced Duffing Oscillator";
% paramNote = "One Oscillator";
% nIC = 10;
% uniformInX = false;
% iAString = "High Attractor"; 

% Two oscillator modications
% rhsString = 'Duffing';
% dim = 4; 
% note = "Two Forced Duffing Oscillator";
% paramNote = "Two Oscillators with 10 nICs";
% nIC = 10;
% uniformInX = true;
% iAString = "High Attractor"; 

% % Three oscillator modications
% rhsString = 'ThreeDuffing';
% dim = 6; 
% note = "Three Forced Duffing Oscillator";
% paramNote = "ThreeDuffing";
% nIC = 10;
% uniformInX = true;
% xcoordinates = false;
% rA = .5; 
% rA1 = .5; 
% iAString = "LHL Attractor"; 

% % N oscillator modications
% rhsString = 'NDuffing';
% dim = 10; 
% note = "Five Forced Duffing Oscillator";
% paramNote = "FiveDuffing";
% nIC = 15;
% uniformInX = true;
% xcoordinates = false;
% rA = .1; 
% rA1 = .2; 


% Test modifications
% methodTest = true;
% saveMemory = true;
% dim = 2;
% note = "methodTest";
% paramNote = "methodTest";
% nIC = 10;
% xcoordinates = false;
