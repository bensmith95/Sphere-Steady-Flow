%% Newton iterative solver
% Solves system until correction < 1e-5 
% Inputs:
% U - guess for U component        V - guess for V component
% W - guess for W component        r - radius vector
% beta - angle/height vector       i - current position index
% alpha - under-relaxation parameter
% Outputs:
% U - solution for U component        
% V - solution for V component
% W - solution for W component
function [U,V,W] = newton(U,V,W,r,beta,i,alpha)
    %% Initialise
    % Space marching parameters
    dbeta = beta(2)-beta(1); Nbeta = length(beta); 
    h1 = r(i)-r(i-1); h2 = r(i)-r(i-2);

    % set BCs
    VBC = V(1,i);

    % Set intitial guess
    Uo = U(:,i);  Vo = V(:,i);  Wo = W(:,i);

    %% Iterate
    J = zeros(3*Nbeta); f = zeros(3*Nbeta,1); q = 1;
    while max(abs(q))>1e-5
        % FINITE DIFFERENCES (excluding boundaries)
        % 2nd order centered differences
        dU_dbeta = (Uo(3:end)-Uo(1:end-2))/(2*dbeta);
        dV_dbeta = (Vo(3:end)-Vo(1:end-2))/(2*dbeta);
        dW_dbeta = (Wo(3:end)-Wo(1:end-2))/(2*dbeta);

        % 1st order backwards differences
        dV_dr = (h2^2*V(2:end-1,i-1)-h1^2*V(2:end-1,i-2)+(h1^2-h2^2)*Vo(2:end-1))/(h1^2*h2-h1*h2^2);
        dW_dr = (h2^2*W(2:end-1,i-1)-h1^2*W(2:end-1,i-2)+(h1^2-h2^2)*Wo(2:end-1))/(h1^2*h2-h1*h2^2);

        % 2nd order centered differences
        d2V_dbeta2 = (Vo(3:end)-2*Vo(2:end-1)+Vo(1:end-2))/dbeta^2;
        d2W_dbeta2 = (Wo(3:end)-2*Wo(2:end-1)+Wo(1:end-2))/dbeta^2;

        % checking how well we satisfy the equations at the given i-th
        % location in the interior points 
        RHS_cont = -(r(i)*dW_dr + 2*Wo(2:end-1) + dU_dbeta);
        RHS_V = d2V_dbeta2 - r(i)^2*Wo(2:end-1).*dV_dr - r(i)*Uo(2:end-1).*dV_dbeta - r(i)*Vo(2:end-1).*Wo(2:end-1);
        RHS_W = d2W_dbeta2 - r(i)^2*Wo(2:end-1).*dW_dr - r(i)*Uo(2:end-1).*dW_dbeta + r(i)*Vo(2:end-1).^2;

        % JACOBIAN CONSTRUCTION
        % beta->inf BCs
        J(1,1+Nbeta*0+[2 1 0]) = [-1 4 -3];
        f(1,1) = -(-3*Uo(1)+4*Uo(2)-Uo(3));

        J(2,1+Nbeta*1) = 1;
        f(2,1) = -Vo(1)+VBC;

        J(3,1+Nbeta*2) = 1;
        f(3,1) = -Wo(1);

        for j=2:Nbeta-1
            % continuity
            J(1+(j-1)*3,j+1+Nbeta*0) = 1/(2*dbeta);
            J(1+(j-1)*3,j-1+Nbeta*0) = -1/(2*dbeta);
            J(1+(j-1)*3,j+Nbeta*2) = r(i)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + 2;
            f(1+(j-1)*3,1) = RHS_cont(j-1);
            % W momentum
            J(2+(j-1)*3,j+Nbeta*0) = r(i)*dW_dbeta(j-1);
            J(2+(j-1)*3,j+Nbeta*1) = -2*r(i)*Vo(j);
            J(2+(j-1)*3,j+Nbeta*2) = r(i)^2*(dW_dr(j-1)+Wo(j)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2)) + 2/dbeta^2;
            J(2+(j-1)*3,j+1+Nbeta*2) = r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            J(2+(j-1)*3,j-1+Nbeta*2) = -r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            f(2+(j-1)*3,1) = RHS_W(j-1);
            % V momentum
            J(3+(j-1)*3,j+Nbeta*0) = r(i)*dV_dbeta(j-1);
            J(3+(j-1)*3,j+Nbeta*1) = r(i)^2*Wo(j)*(h1^2-h2^2)/(h1^2*h2-h1*h2^2) + r(i)*Wo(j) + 2/dbeta^2;
            J(3+(j-1)*3,j+1+Nbeta*1) = r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            J(3+(j-1)*3,j-1+Nbeta*1) = -r(i)*Uo(j)/(2*dbeta) - 1/dbeta^2;
            J(3+(j-1)*3,j+Nbeta*2) = r(i)^2*dV_dr(j-1) + r(i)*Vo(j);
            f(3+(j-1)*3,1) = RHS_V(j-1); 
        end

        % beta=0 BCs
        J(end-2,Nbeta*1) = 1;
        f(end-2,1) = -Uo(end);

        J(end-1,Nbeta*2-[2 1 0]) = [1 -4 3];
        f(end-1,1) = -(3*Vo(end)-4*Vo(end-1)+Vo(end-2));

        J(end,Nbeta*3-[2 1 0]) = [1 -4 3];
        f(end,1) = -(3*Wo(end)-4*Wo(end-1)+Wo(end-2));

        % solve for correction, q
        q = J\f;

        %fprintf('Max q = %.6f\n',max(abs(q))); 
        % Newton-Raphson iteration
        Uo = Uo+q([1:Nbeta]+Nbeta*0)*alpha; 
        Vo = Vo+q([1:Nbeta]+Nbeta*1)*alpha;
        Wo = Wo+q([1:Nbeta]+Nbeta*2)*alpha; 
    end     
    % set solution
    U(:,i) = Uo; V(:,i) = Vo; W(:,i) = Wo; 
end