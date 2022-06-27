function [J] = DuffingFilteredJacobian(qo,Mrhs)
        A = [0 1 0;...
            -3*Mrhs.a3*qo(1)^2-Mrhs.a1, -Mrhs.nu 1;
            0 0 -Mrhs.wc]; %linearized original system at attractor coordinate
        E = [0 0 0;0 0 0; 0 0 Mrhs.wc^2]; %contribution of p to original system
        J = [A E;zeros(size(E)) -A.']; %Jacobian of hamiltonian system 
end

