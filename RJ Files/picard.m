%% Picard iterative solver
% Solves system until correction < 0.1 
% in order to provide a close enough guess for Newton method
% Inputs:
% U - guess for U component        V - guess for V component
% W - guess for W component        r - radius vector
% beta - angle/height vector       i - current position index
% alpha - under-relaxation parameter
% Outputs:
% U - improved guess for U component        
% V - improved guess for V component
% W - improved guess for W component

function [U,V,W] = picard(U,V,W,r,beta,i,alpha)
    %% initialise
    % Space marching parameters
    dbeta = beta(2)-beta(1); Nbeta = length(beta); 
    h1 = r(i)-r(i-1); h2 = r(i)-r(i-2);

    % set BCs
    VBC = V(1,i);

    % Set intitial guess
    Uo = U(:,i-1);  Vo = V(:,i-1);  Wo = W(:,i-1);
    
    %% iterate
    A = zeros(3*Nbeta); b = zeros(3*Nbeta,1); q = 1;
    while max(abs(q))>1e-1

        % beta->Inf BCs
        A(1,1+Nbeta*0+[2 1 0]) = [-1 4 -3]; % U
        A(2,1+Nbeta*1) = 1; b(2) = VBC; % V
        A(3,1+Nbeta*2) = 1; % W

        for j=2:Nbeta-1
            % Continuity
            A(1+(j-1)*3,j-1+Nbeta*0) = -1/(2*dbeta);
            A(1+(j-1)*3,j+1+Nbeta*0) = 1/(2*dbeta);
            A(1+(j-1)*3,j+Nbeta*2) = r(i)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + 2; 
            b(1+(j-1)*3) = r(i)*(h1^2*W(j,i-2)-h2^2*W(j,i-1))/(h1^2*h2-h1*h2^2);
            % V Momentum
            A(2+(j-1)*3,j-1+Nbeta*1) = -r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            A(2+(j-1)*3,j+Nbeta*1) = r(i)^2*Wo(j)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + r(i)*Wo(j) + 2/dbeta^2;
            A(2+(j-1)*3,j+1+Nbeta*1) = r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            b(2+(j-1)*3) = r(i)^2*Wo(j)*(h1^2*V(j,i-2)-h2^2*V(j,i-1))/(h1^2*h2-h1*h2^2);
            % W Momentum
            A(3+(j-1)*3,j+Nbeta*1) = -r(i)*Vo(j);
            A(3+(j-1)*3,j-1+Nbeta*2) = -r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            A(3+(j-1)*3,j+Nbeta*2) = r(i)^2*Wo(j)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + 2/dbeta^2;
            A(3+(j-1)*3,j+1+Nbeta*2) = r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            b(3+(j-1)*3) = r(i)^2*Wo(j)*(h1^2*W(j,i-2)-h2^2*W(j,i-1))/(h1^2*h2-h1*h2^2);
        end

        % beta=0 BCs
        A(end-2,Nbeta*1) = 1; % U
        A(end-1,Nbeta*2-[2 1 0]) = [1 -4 3]; % V
        A(end,Nbeta*3-[2 1 0]) = [1 -4 3]; % W

        % Solve system
        s = A\b;

        % Under-relax
        Us = (1-alpha)*Uo+s([1:Nbeta]+Nbeta*0)*alpha;
        Vs = (1-alpha)*Vo+s([1:Nbeta]+Nbeta*1)*alpha;
        Ws = (1-alpha)*Wo+s([1:Nbeta]+Nbeta*2)*alpha;

        % Compute change in solution
        Q = abs([Us-Uo;Vs-Vo;Ws-Wo]);
        q = max(Q,[],'all');

        % set new guess
        Uo = Us; Vo = Vs; Wo = Ws;

    end
    % set new solution
    U(:,i) = Uo; V(:,i) = Vo; W(:,i) = Wo;
end