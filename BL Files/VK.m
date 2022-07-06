%% Solves Von Karman flow for rotating disk
% Inputs:
% Zmax - maximum value to solve along
% N - number of discretisation points
% Outputs:
% Vel - cell array of velocities and pressure distribution
% z - vector which von karman equations are solved upon
function [Vel,z] = VK(Zmax,N)
    % determine log spacing, so better precision near sphere surface
    z = [0,logspace(-2,log10(Zmax),N)];

    % solve von karman system
    options = bvpset('RelTol',1e-8,'AbsTol',1e-10,'NMax',10^4); % set tolerances
    solinit = bvpinit(z,[1 0 1 0 1]); % set initial guess
    % numerical solver
    Y = bvp4c(@(x,U) VK_system(x,U),@(Ua,Ub) VK_BC(Ua,Ub),solinit,options);
    Sol = deval(Y,z)';

    % save solution to cell array
    for i = 1:5
        Vel{i} = Sol(:,i);
    end
end

%% Von Karman/O(1) Banks system
function dudz = VK_system(x,U)
    dudz = [U(2); U(1)^2-U(3)^2+U(5)*U(2);
            U(4); U(5)*U(4)+2*U(1)*U(3);
            -2*U(1)];
end
%% Set boundary conditions
function BC = VK_BC(Ua,Ub)
    BC = [Ua(1);Ub(1);Ua(3)-1;Ub(3);Ua(5)];
end