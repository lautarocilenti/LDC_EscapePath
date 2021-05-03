function [M] = Parameters(parameterNames,parameterValues)
%Configure PARAMETERS here 
note = "Forced Duffing Oscillator";
paramNote = "Duffing";
a1 = 1; a3 = .3; nu = .1; F = .4; w = 1.4; %rhs parameters (note basin interpolant mat file must be changed if rhs parameters are changed)
rIC = 10^-6; %radius of momenta initial conditions
qo = [1.3590;2.4170]; %initial condition in phase space
% qo = [-sqrt(-a1/a3);0]; %initial condition in phase space for bistable
BN = 1; %Initial Basin Identifier
nIC = 5; %number of initial conditions in circle around qo
rhsString = 'Duffing';
 T = 2*pi/w;  dt = T/32; tf = 300*T;
solver = @ode45;
psiEps = .05; %phase threshold


%Calculated parameters
D = 2*length(qo);
tspan = [0:dt:tf];
dtheta =(2*pi)/(nIC);
theta = [dtheta:dtheta:2*pi];

%higher resolution
ub = 1.2; lb = .9;
dtheta2 =(ub-lb)/(nIC);
theta2 = lb+dtheta2:dtheta2:ub;
nIC = 2*nIC;
theta = [theta,theta2];
theta = sort(theta,'ascend')

% nIC = 1;
% theta = 1.0598;



%Store Parameters in a structure
M.note = note;
M.rIC = rIC; M.qo = qo; M.BN  = BN;
M.nIC = nIC;
M.rhsString = rhsString;
M.tf = tf; M.dt = dt;
M.solver = solver;
M.D = D;
M.tspan = tspan;
M.theta = theta; M.dtheta = dtheta;

M.paramNote = paramNote;
M.Mrhs.a1 = a1; M.Mrhs.a3 = a3; M.Mrhs.nu = nu;
M.Mrhs.F = F; M.Mrhs.w = w; M.Mrhs.psiEps = psiEps;


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

