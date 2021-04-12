function [M] = Parameters(parameterNames,parameterValues)
%Configure PARAMETERS here 
note = "Overdamped, Bistable Unforced Duffing Oscillator";
paramNote = "OverDamped";
a1 = -1; a3 = 1; nu = 1; %rhs parameters (note basin interpolant mat file must be changed if rhs parameters are changed)
rIC = 10^-6; %radius of momenta initial conditions
qo = [-sqrt(-a1/a3);0]; %initial condition in phase space
BN = 1; %Initial Basin Identifier
nIC = 100; %number of initial conditions in circle around qo
rhsString = 'Unforced';
tf = 400; dt = 1E-2; 
solver = @ode45;
gsTrials = 1000; %total trials
gsInTrials = 200; %initial trials


%Calculated parameters
D = 2*length(qo);
tspan = [0,tf];
dtheta =(2*pi)/(nIC);
theta = [dtheta:dtheta:2*pi];

%higher resolution
ub = 1.9; lb = 1.83;
dtheta2 =(ub-lb)/(nIC);
theta2 = lb+dtheta2:dtheta2:ub;
nIC = 2*nIC;
theta = [theta,theta2];
theta = sort(theta,'ascend')



%Store Parameters in a structure
M.note = note;
M.rIC = rIC; M.qo = qo; M.BN  = 1;
M.nIC = nIC;
M.rhsString = rhsString;
M.tf = tf; M.dt = dt;
M.solver = solver;
M.D = D;
M.tspan = tspan;
M.theta = theta; M.dtheta = dtheta;
M.gsTrials = gsTrials; M . gsInTrials = gsInTrials;

M.paramNote = paramNote;
M.Mrhs.a1 = a1; M.Mrhs.a3 = a3; M.Mrhs.nu = nu;


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


%Set Moment RHS
M.RHS = str2func([M.rhsString,'RHS']);
M.Lagrangian = str2func([M.rhsString,'Lagrangian'])


end

