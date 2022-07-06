%% Solves the boundary layer equations 
function BLRegion
%% Pre-Condition
    % Space marching parameters
    eta_max = 30; Neta = 301; Ntheta = 1501;
    eta = linspace(0,eta_max,Neta); deta = eta(2)-eta(1);
    theta = linspace(0,pi/2,Ntheta); dtheta = theta(2)-theta(1);

    % Solve Von Karman equations
    str1 = fprintf('Solving Von Karman equations...\n');
    [Vel,z] = VK(eta_max+1,Neta);
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solved Von Karman equations.\n');

    % assign velocity profiles
    U1 = Vel{1}; V1 = Vel{3}; W1 = Vel{5}; 

    % initialise field
    % field discretised as theta x eta, or M_ij = [theta(i),eta(j)]
    U = repmat(theta',1,Neta).*repmat(spline(z,U1,eta),Ntheta,1);
    V = repmat(theta',1,Neta).*repmat(spline(z,V1,eta),Ntheta,1);
    W = repmat(spline(z,W1,eta),Ntheta,1);

%% NEWTON RAPHSON MULTIVARIATE SCHEME
    str2 = fprintf('Solving Boundary Layer equations at theta =  %.2f\n',0);
    for i=3:Ntheta
        deg = theta(i)*180/pi; fprintf(repmat('\b',1,str2)); 
        str2 = fprintf('Solving Boundary Layer equations at theta =  %.2f\n',deg);

        % Use previous latitude as an intitial guess
        Uo = U(i-1,:);  Vo = V(i-1,:);  Wo = W(i-1,:);

        % Iterate
        gg = 0.95; q = 1;
        A = zeros(3*Neta); B = zeros(3*Neta,1);
        while max(abs(q))>1e-5
            % FINITE DIFFERENCES (excluding boundaries)
            % 2nd order centered differences
            dU_deta = (Uo(3:end)-Uo(1:end-2))/(2*deta);
            dV_deta = (Vo(3:end)-Vo(1:end-2))/(2*deta);
            dW_deta = (Wo(3:end)-Wo(1:end-2))/(2*deta);            

            % 2nd order three point backwards differences
            dU_dtheta = (3*Uo(2:end-1)-4*U(i-1,2:end-1)+U(i-2,2:end-1))/(2*dtheta);
            dV_dtheta = (3*Vo(2:end-1)-4*V(i-1,2:end-1)+V(i-2,2:end-1))/(2*dtheta);

            % 2nd order centered differences
            d2V_deta2 = (Vo(3:end)-2*Vo(2:end-1)+Vo(1:end-2))/deta^2;
            d2U_deta2 = (Uo(3:end)-2*Uo(2:end-1)+Uo(1:end-2))/deta^2;

            % checking how well we satisfy the equations at the given i-th
            % location in the interior points
            RHS_cont = -(dU_dtheta + dW_deta + Uo(2:end-1)*cot(theta(i)));
            RHS_U = d2U_deta2 - Uo(2:end-1).*dU_dtheta - Wo(2:end-1).*dU_deta + Vo(2:end-1).^2*cot(theta(i));
            RHS_V = d2V_deta2 - Uo(2:end-1).*dV_dtheta - Wo(2:end-1).*dV_deta - Uo(2:end-1).*Vo(2:end-1)*cot(theta(i));

            % JACOBIAN CONSTRUCTION
            % wall boundary conditions
            A(1,1+Neta*0) = 1;
            B(1,1) = -Uo(1);

            A(2,1+Neta*1) = 1;
            B(2,1) = -Vo(1) + sin(theta(i));

            A(3,1+Neta*2) = 1;
            B(3,1) = -Wo(1);

            % internal points
            for j=2:Neta-1
                % continuity
                A(1+(j-1)*3,j+Neta*0) = 3/(2*dtheta) + cot(theta(i));
                A(1+(j-1)*3,j+1+Neta*2) = 1/(2*deta);
                A(1+(j-1)*3,j-1+Neta*2) = -1/(2*deta);
                B(1+(j-1)*3,1) = RHS_cont(j-1);
                % U momentum
                A(2+(j-1)*3,j+Neta*0) = dU_dtheta(j-1) + 3*Uo(j)/(2*dtheta) - (-2/deta^2);
                A(2+(j-1)*3,j+1+Neta*0) = Wo(j)/(2*deta) - 1/deta^2;
                A(2+(j-1)*3,j-1+Neta*0) = -Wo(j)/(2*deta) - 1/deta^2;
                A(2+(j-1)*3,j+Neta*1) = -2*Vo(j)*cot(theta(i));
                A(2+(j-1)*3,j+Neta*2) = dU_deta(j-1);
                B(2+(j-1)*3,1) = RHS_U(j-1);
                % V momentum
                A(3+(j-1)*3,j+Neta*0) = dV_dtheta(j-1) + Vo(j)*cot(theta(i));
                A(3+(j-1)*3,j+Neta*1) = 3*Uo(j)/(2*dtheta) + Uo(j)*cot(theta(i)) - (-2/deta^2);
                A(3+(j-1)*3,j+1+Neta*1) = Wo(j)/(2*deta) - 1/deta^2;
                A(3+(j-1)*3,j-1+Neta*1) = -Wo(j)/(2*deta) - 1/deta^2;
                A(3+(j-1)*3,j+Neta*2) = dV_deta(j-1);       
                B(3+(j-1)*3,1) = RHS_V(j-1);  
            end

            % top boundary condition        
            A(end-2,Neta) = 1;
            B(end-2,1) = -Uo(end);

            A(end-1,Neta*2) = 1;
            B(end-1,1) = -Vo(end);

            A(end,3*Neta-[2 1 0]) = [1 -4 3];
            B(end,1) = -(3*Wo(end)-4*Wo(end-1)+Wo(end-2));

            % solve for correction, q
            q = A\B;
            %fprintf('Max q = %.6f\n',max(q))
            % Newton-Raphson iteration
            Uo = Uo+q([1:Neta]+Neta*0)'*gg;
            Vo = Vo+q([1:Neta]+Neta*1)'*gg;
            Wo = Wo+q([1:Neta]+Neta*2)'*gg;        
        end

        % set solution
        U(i,:) = Uo; V(i,:) = Vo; W(i,:) = Wo;
    end
    fprintf(repmat('\b',1,str2)); str2 = fprintf('All latitudes solved for.\n');
    
%% Post Process 
    str3 = fprintf('Post-Processing...\n');

    % Determine pressure
    dP = U.^2+V.^2;
    P = cumtrapz(flip(eta),flip(dP,2),2); P = flip(P,2);

    % determine 2nd order finite differences for derivatives 
    dU_dtheta = [-(3*U(1,:)-4*U(2,:)+U(3,:)); U(3:end,:)-U(1:end-2,:); (3*U(end,:)-4*U(end-1,:)+U(end-2,:))]/2/dtheta;
    dV_dtheta = [-(3*V(1,:)-4*V(2,:)+V(3,:)); V(3:end,:)-V(1:end-2,:); (3*V(end,:)-4*V(end-1,:)+V(end-2,:))]/2/dtheta;
    dW_dtheta = [-(3*W(1,:)-4*W(2,:)+W(3,:)); W(3:end,:)-W(1:end-2,:); (3*W(end,:)-4*W(end-1,:)+W(end-2,:))]/2/dtheta;
    dP_dtheta = [-(3*P(1,:)-4*P(2,:)+P(3,:)); P(3:end,:)-P(1:end-2,:); (3*P(end,:)-4*P(end-1,:)+P(end-2,:))]/2/dtheta;

    dU = [-(3*U(:,1)-4*U(:,2)+U(:,3)), U(:,3:end)-U(:,1:end-2), (3*U(:,end)-4*U(:,end-1)+U(:,end-2))]/2/deta;
    dV = [-(3*V(:,1)-4*V(:,2)+V(:,3)), V(:,3:end)-V(:,1:end-2), (3*V(:,end)-4*V(:,end-1)+V(:,end-2))]/2/deta;
    dW = [-(3*W(:,1)-4*W(:,2)+W(:,3)), W(:,3:end)-W(:,1:end-2), (3*W(:,end)-4*W(:,end-1)+W(:,end-2))]/2/deta;

    % organise data into a cell array
    clear Vel
    VelBL{1} = U; VelBL{2} = V; VelBL{3} = W; VelBL{4} = P; 
    VelBL{5} = dU; VelBL{6} = dV; VelBL{7} = dW; VelBL{8} = dP;  
    VelBL{9} = dU_dtheta; VelBL{10} = dV_dtheta; VelBL{11} = dW_dtheta; VelBL{12} = dP_dtheta; 

    % save data to file
    filename = '../Flows/BL.mat';
    if exist('../Flows','dir')==0; mkdir ../Flows; end
    save(filename, 'VelBL', 'eta', 'theta')
    fprintf(repmat('\b',1,str3)); str3 = fprintf('Flow saved in %s\n', filename); pause(1)
    fprintf(repmat('\b',1,str3));fprintf(repmat('\b',1,str2));fprintf(repmat('\b',1,str1));
end