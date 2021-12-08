function [xoSet,M] = GenerateInitialConditions(thetaSet,M)
%GENERATEINITIALCONDITIONS Creates set of initial conditions
%
% Input: thetaSet: set of angles surrounding fixed point of poincare section
%        M: large  structure of parameters
%
% Output: xoSet: each column contains 2D values and is an initial condition
%                of the 2D hamiltonian system
% 


if isempty(M.transMatrix)
    %objects for rhs generation 
        x = sym('x',[M.dim*2,1]);
        syms tau

    %Grab the location of the attractor in the poincare section
    [qo] = GetInitialFixedPoint(M.pp,M.Mrhs);


    if strcmp(M.rhsString,'Duffing') %Duffing forced case
        A = [0 1; -3*M.Mrhs.a3*qo(1)^2-M.Mrhs.a1, -M.Mrhs.nu]; %linearized original system at attractor coordinate
        E = [0 0;0 1]; %contribution of p to original system
        J = [A E;zeros(size(E)) -A.']; %Jacobian of hamiltonian system 
    elseif strcmp(M.rhsString,'TwoDuffing')
        J = TwoDuffingJacobian(qo(1),qo(3),M.Mrhs); 
    elseif strcmp(M.rhsString,'ThreeDuffing')
        J = ThreeDuffingJacobian([qo(1),qo(3),qo(5)],M.Mrhs); 
    elseif strcmp(M.rhsString,'NDuffing')
        Xj = [qo',zeros(size(qo'))];
        J = NDuffingJacobian(Xj,M.Mrhs); 
    else
        error("Unknown Jacobian\n")
    end

    jRHS = matlabFunction(J*x,'vars',{tau,x}); %convert to function handle for integration

    tspan = [0,M.Mrhs.T]; %one period time span
    I = eye(M.dim*2); %identity matrix
    opts = odeset('RelTol',1e-12,'AbsTol',1e-12); %ode45 options
    Psi_T = zeros(size(I)); %initialize Psi

    for i = 1:M.dim*2
        xo = I(:,i); %initial condition is identity matrix columns
        [t,y] =  M.solver(jRHS, tspan, xo,[opts]);
        Psi_T(:,i) = y(end,:).'; %columns of Psi at time T
    end

    E = Psi_T;

    %get eigenvalues
    [ev,e] = eig(E);

    lambda = diag(e); %eigenvalues 
    [lambda,isortlambda] = sort(vecnorm(lambda,2,2),'ascend'); %sort eigenvalues by distance from origin in real imag plane
    iUnstableLambda = find(lambda>1); 
    iStableLambda = find(lambda<=1);

    ev = ev(:,isortlambda); %re sort eigenvectors 

    if M.dim == 2
        u_ev = ev(:,iUnstableLambda(1)); %eigenvectors of unstable eigenvalues, just grab 1 of conjugate pair
    elseif M.dim ==4
        u_ev = ev(:,iUnstableLambda([1,3])); %eigenvectors of unstable eigenvalues, just grab 1 per conjugate pair
    elseif M.dim == 6
        u_ev = ev(:,iUnstableLambda([1,3,5])); %eigenvectors of unstable eigenvalues, just grab 1 per conjugate pair
    elseif M.dim > 6
        iev = 1:2:M.dim;
        u_ev = ev(:,iUnstableLambda([iev])); %eigenvectors of unstable eigenvalues, just grab 1 per conjugate pair
    end
    Z_xv = [real(u_ev(1:M.dim,:)),imag(u_ev(1:M.dim,:))];
    Z_pxpv = [real(u_ev(M.dim+1:2*M.dim,:)),imag(u_ev(M.dim+1:2*M.dim,:))];
    T = Z_xv\Z_pxpv; %transformation matrix
    M.transMatrix = T;
    
end
%     

if M.xcoordinates
        q_eps = M.rIC*thetaSet;
else

    if M.dim == 2
        q_eps = [M.rIC*cos(thetaSet);M.rIC*sin(thetaSet)]; %initial offset from attractor
    elseif M.dim ==4
        if contains(M.searchAlgorithm,"OneO")
            if M.io == 1
                q_eps = [M.rIC*cos(thetaSet);zeros(size(thetaSet));M.rIC*sin(thetaSet);zeros(size(thetaSet))]; 
            else
                q_eps = [zeros(size(thetaSet));M.rIC*cos(thetaSet);zeros(size(thetaSet));M.rIC*sin(thetaSet)]; 
            end
            
        else
            theta1 = thetaSet(1,:);
            theta2 = thetaSet(2,:);
            theta3 = thetaSet(3,:);
            q_eps = [M.rIC*cos(theta1);M.rIC*sin(theta1).*cos(theta2);M.rIC*sin(theta1).*sin(theta2).*cos(theta3);M.rIC*sin(theta1).*sin(theta2).*sin(theta3)]; %initial offset from attractor 
        end
    elseif M.dim == 6
        q_eps = [M.rIC*cos(thetaSet(1,:)); ...
                 M.rIC*sin(thetaSet(1,:)).*cos(thetaSet(2,:)); ...
                 M.rIC*sin(thetaSet(1,:)).*sin(thetaSet(2,:)).*cos(thetaSet(3,:)); ... 
                 M.rIC*sin(thetaSet(1,:)).*sin(thetaSet(2,:)).*sin(thetaSet(3,:)).*cos(thetaSet(4,:)); ... 
                 M.rIC*sin(thetaSet(1,:)).*sin(thetaSet(2,:)).*sin(thetaSet(3,:)).*sin(thetaSet(4,:)).*cos(thetaSet(5,:)); ... 
                 M.rIC*sin(thetaSet(1,:)).*sin(thetaSet(2,:)).*sin(thetaSet(3,:)).*sin(thetaSet(4,:)).*sin(thetaSet(5,:))]; %initial offset from attractor 
     elseif M.dim > 6
         % build n sphere from angles
         q_eps = zeros(M.dim,size(thetaSet,2));
         for i = 1:M.dim
             if i == M.dim
                jSines = 1:i-1;
                jCos = [];
             else
                jSines = 1:i-1;
                jCos = i;
             end
            if ~isempty(jSines)
                sines = sin(thetaSet(jSines,:));
                sines = prod(sines,1);
            else
                sines = ones(1,size(thetaSet,2));
            end
            if ~isempty(jCos)
            	cosines = cos(thetaSet(jCos,:));
            else
                 cosines = ones(1,size(thetaSet,2));
            end
            q_eps(i,:) = M.rIC.*sines.*cosines;
         end
    end
end


T = M.transMatrix;
p_eps = T*q_eps; %initial momenta

eps = ArrangeCoordinates([q_eps;p_eps],M);

xoSet = [qo.*ones(M.dim,size(thetaSet,2));zeros(M.dim,size(thetaSet,2))]+eps;
    


% save('temp.mat','e','ev','q_eps')



end

function [qo] = GetInitialFixedPoint(pp,Mrhs)
    t = pp/Mrhs.w;
    if t~=0
        error("need to interpolate fixed point")
    end
    [fp] = Mrhs.FixedPoints.FP(Mrhs.iA,:);
    qo = fp';
end

function [y] = ArrangeCoordinates(x,M)
    y = x;

    if M.dim == 4
        % y = [x1 v1 x2 v2 p1 p1d p2 p2d]; %hamiltonian rhs
        % x = [x1 x2 v1 v2 p1 p2 p1d p2d]; %jacobian rhs
        y(2,:) = x(3,:);
        y(3,:) = x(2,:);
        y(6,:) = x(7,:);
        y(7,:) = x(6,:);
    elseif M.dim == 6
        % y = [x1 v1 x2 v2 x3 v3 p1 p1d p2 p2d p3 p3d]; %hamiltonian rhs
        % x = [x1 x2 x3 v1 v2 v3 p1 p2 p3 p1d p2d p3d]; %jacobian rhs
        y(2,:) = x(4,:);
        y(3,:) = x(2,:);
        y(4,:) = x(5,:);
        y(5,:) = x(3,:);
        y(8,:) = x(10,:);
        y(9,:) = x(8,:);
        y(10,:) = x(11,:);
        y(11,:) = x(9,:);  
    elseif M.dim >6
        %x =[ix iv ipx ipv]
        L = size(x,1)/2;
        ix1 = 1:L/2;
        iv1 = L/2+1:L;
        ipx1 = L+1:L+L/2;
        ipv1 = L+L/2+1:2*L;
        ix2 = 1:2:L;
        iv2 = 2:2:L;
        ipx2 = L+1:2:2*L;
        ipv2 = L+2:2:2*L;
        y(ix2) = x(ix1);
        y(iv2) = x(iv1);
        y(ipx2) = x(ipx1);
        y(ipv2) = x(ipv1);
        
    end

end