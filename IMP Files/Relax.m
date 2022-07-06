%% Gauss-Seidel-Newton block iterative scheme to relax/solve SF-Vort. equations
% Inputs:
% Psi - Approximation for SF           W - Approximation for Vort.
% V - Approximation for azi. vel.      Sp - RHS of SF eqn                  
% Sw - RHS of Vort. eqn                Sv - RHS of V eqn
% BC - Bound. Cond. function           h - size of grid
% Re - Reynolds Number                 n - # of iterations
% a - under-relaxation parameter
% Outputs:
% Psi - Updated approx. for SF        
% W - Updated approx. for Vort.
% V - Updated approx for azi. vel.     
function [Psi,W,V] = Relax(Psi,W,V,Sp,Sw,Sv,BC,h,Re,n,a)
    %% Initalise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2; 
    eps = 1/sqrt(Re); hi = 1/h; [Pxxt,~] = BC(h,Re);
    
    %% Iterate
    % direction of iterations follows overall direction of flow
    for k = 1:n
        for j = Ny:-1:2
            for i = 3:Nx+1
                dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi; 
                dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;
                dVdy = (V(j+1,i)-V(j-1,i))*0.5*hi; 
                if j>2
                    % Internal Points
                    if i<Nx+1
                        A = zeros(3); b = zeros(3,1);

                        % SF & V mom. upwind terms
                        if dPsidy>=0
                            dWdx = (3*W(j,i)-4*W(j,i-1)+W(j,i-2))*0.5*hi; 
                            dV_dx = (3*V(j,i)-4*V(j,i-1)+V(j,i-2))*0.5*hi;
                            X = 3*dPsidy*0.5*hi;
                        else
                            dWdx = -(3*W(j,i)-4*W(j,i+1)+W(j,i+2))*0.5*hi; 
                            dV_dx = -(3*V(j,i)-4*V(j,i+1)+V(j,i+2))*0.5*hi;
                            X = -3*dPsidy*0.5*hi;
                        end
                        if -dPsidx>=0
                            dWdy = (3*W(j,i)-4*W(j-1,i)+W(j-2,i))*0.5*hi; 
                            dV_dy = (3*V(j,i)-4*V(j-1,i)+V(j-2,i))*0.5*hi;
                            Y = -3*dPsidx*0.5*hi;
                        else
                            dWdy = -(3*W(j,i)-4*W(j+1,i)+W(j+2,i))*0.5*hi; 
                            dV_dy = -(3*V(j,i)-4*V(j+1,i)+V(j+2,i))*0.5*hi;
                            Y = 3*dPsidx*0.5*hi;
                        end
                        LapW = (W(j-1,i)+W(j,i-1)-4*W(j,i)+W(j,i+1)+W(j+1,i))*hi^2;
                        LapV = (V(j-1,i)+V(j,i-1)-4*V(j,i)+V(j,i+1)+V(j+1,i))*hi^2;

                        % SF
                        A(1,2) = X+Y+4*eps*hi^2; A(1,3) = 2*dVdy;
                        b(1) = -(dPsidy*dWdx - dPsidx*dWdy + 2*V(j,i)*dVdy - eps*LapW - Sp(j,i));

                        % Vorticity
                        LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;
                        A(2,1) = -4*hi^2; A(2,2) = 1;
                        b(2) = -(W(j,i) + LapPsi - Sw(j,i));

                        % V momentum
                        A(3,3) = X+Y+4*eps*hi^2;
                        b(3) = -(dPsidy*dV_dx - dPsidx*dV_dy - eps*LapV - Sv(j,i));

                        % solve & correct
                        s = A\b;
                        Psi(j,i) = Psi(j,i) + a*s(1); 
                        W(j,i) = W(j,i) + a*s(2);
                        V(j,i) = V(j,i) + a*s(3);

                        % Update BC's 
                        if ismember(i,[3,4])
                            W(j,2) = -2*Psi(j,3)*hi^2; W(j,1) = 3*W(j,2) - 3*W(j,3) + W(j,4); 
                            Psi(j,1) = Psi(j,3); 
                            V(j,1) = 3*V(j,2) - 3*V(j,3) + V(j,4);
                        end
                        if ismember(i,[Nx-1,Nx])
                            if j<Ny; W(j+1,Nx+1) = (2*Psi(j+1,Nx+1)-Psi(j,Nx)-Psi(j+2,Nx))*hi^2; end
                            if j>3; W(j-1,Nx+1) = (2*Psi(j-1,Nx+1)-Psi(j-2,Nx)-Psi(j,Nx))*hi^2; end
                            W(j,Nx+2) = 3*W(j,Nx+1) - 3*W(j,Nx) + W(j,Nx-1); Psi(j,Nx+2) = Psi(j,Nx); 
                            V(j,Nx+2) = V(j,Nx); 
                        end
                        if ismember(j,[3,4]) 
                            W(1,i) = -3*W(3,i) + W(4,i); Psi(1,i) = 3*Psi(2,i) - 3*Psi(3,i) + Psi(4,i); 
                            V(1,i) = V(3,i);
                        end
                        if ismember(j,[Ny-1,Ny])
                            W(Ny+1,i) = 2*(Psi(Ny+1,i)-Psi(Ny,i))*hi^2 - Pxxt(i-1)'; W(Ny+2,i) = 3*W(Ny+1,i) - 3*W(Ny,i) + W(Ny-1,i); 
                            Psi(Ny+2,i) = Psi(Ny,i); 
                            V(Ny+2,i) = 3*V(Ny+1,i) - 3*V(Ny,i) + V(Ny-1,i);
                        end
                    % Neumann BC
                    else
                        Psi(j,Nx+1) = -(Psi(j-1,Nx+1)-Psi(j-1,Nx)+Psi(j,Nx)+Psi(j,Nx+2)+Psi(j+1,Nx+1)-Psi(j+1,Nx)-Sp(j,Nx+1))/(-2);
                        W(j,Nx+1) = (2*Psi(j,Nx+1)-Psi(j-1,Nx)-Psi(j+1,Nx))*hi^2;
                        V(j,Nx+1) = (V(j-1,Nx+1)+V(j,Nx)+V(j,Nx+2)+V(j+1,Nx+1)-Sv(j,Nx+1))/4;
                                                
                        % BCs
                        W(j,Nx+2) = 3*W(j,Nx+1) - 3*W(j,Nx) + W(j,Nx-1);
                        if j==3; Psi(1,Nx+1) = Psi(3,Nx+1); V(1,Nx+1) = V(3,Nx+1); end 
                    end
                else
                    if i<Nx+1
                        % V mom
                        if dPsidy>=0
                            dV_dx = (3*V(2,i)-4*V(2,i-1)+V(2,i-2))*0.5*hi; 
                            X = 3*dPsidy*0.5*hi;
                        else
                            dV_dx = -(3*V(2,i)-4*V(2,i+1)+V(2,i+2))*0.5*hi; 
                            X = -3*dPsidy*0.5*hi;
                        end
                        LapV = (V(1,i)+V(2,i-1)-4*V(2,i)+V(2,i+1)+V(3,i))*hi^2;  
                        V(2,i) = V(2,i) - a*(dPsidy*dV_dx - eps*LapV - Sv(2,i))/(X + 4*eps*hi^2);
                        
                        % BCs 
                        if ismember(i,[3,4]); V(2,1) = 3 - 3*V(2,3) + V(2,4); end
                        if i==Nx; V(2,Nx+2)=V(2,Nx); end
                    else
                        V(2,Nx+1) = (V(1,Nx+1)+V(2,Nx)+V(2,Nx+2)+V(3,Nx+1)-Sv(2,Nx+1))/4;
                    end
                end
            end
        end
        % Update BC's
        W(Ny+1,2:Nx+1) = 2*(Psi(Ny+1,2:Nx+1)-Psi(Ny,2:Nx+1))*hi^2 - Pxxt(:)';
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(3:Ny,Nx+1) = (2*Psi(3:Ny,Nx+1)-Psi(2:Ny-1,Nx)-Psi(4:Ny+1,Nx))*hi^2;
        % Extrapolated points
        Psi(2:Ny+1,1) = Psi(2:Ny+1,3);
        Psi(2:Ny+1,Nx+2) = Psi(2:Ny+1,Nx);
        Psi(1,2:Nx+1) = 3*Psi(2,2:Nx+1) - 3*Psi(3,2:Nx+1) + Psi(4,2:Nx+1);
        Psi(Ny+2,2:Nx+1) = Psi(Ny,2:Nx+1);
        W(2:Ny+1,1) = 3*W(2:Ny+1,2) - 3*W(2:Ny+1,3) + W(2:Ny+1,4);
        W(2:Ny+1,Nx+2) = 3*W(2:Ny+1,Nx+1) - 3*W(2:Ny+1,Nx) + W(2:Ny+1,Nx-1);
        W(1,2:Nx+1) = -3*W(3,2:Nx+1) + W(4,2:Nx+1);
        W(Ny+2,2:Nx+1) = 3*W(Ny+1,2:Nx+1) - 3*W(Ny,2:Nx+1) + W(Ny-1,2:Nx+1);
        V(2:Ny+1,1) = 3*V(2:Ny+1,2) - 3*V(2:Ny+1,3) + V(2:Ny+1,4);
        V(2:Ny+1,Nx+2) = V(2:Ny+1,Nx);
        V(1,2:Nx+1) = V(3,2:Nx+1);
        V(Ny+2,2:Nx+1) = 3*V(Ny+1,2:Nx+1) - 3*V(Ny,2:Nx+1) + V(Ny-1,2:Nx+1);
    end
end