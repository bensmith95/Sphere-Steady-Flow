%% Determines new Right Hand Side of equations for V-Cycles
% Inputs:
% Psi - Approximation for SF           W - Approximation for Vort.
% V - Approximation for azi. vel.      Pres - residual of SF eqn            
% Wres - residual of Vort. eqn         Vres - residual of V eqn
% h - size of grid                     Re - Reynolds Number 
% Outputs:
% Sp - new RHS of SF equation          
% Sw - new RHS of Vort. eqn
% Sv - new RHS of V eqn                
function [Sph,Swh,Svh] = RHS(Psi,W,V,Pres,Wres,Vres,h,Re)
    %% Initialise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2;
    Sph = zeros(N); Swh = zeros(N); Svh = zeros(N);
    
    %% Determine A(v2h)
    eps = 1/sqrt(Re); hi = 1/h;
    for j = Ny:-1:2
        for i = 3:Nx+1
            dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi; 
            dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;
            dVdy = (V(j+1,i)-V(j-1,i))*0.5*hi; 
            if j>2
                if i<Nx+1
                    % SF & V terms
                    if dPsidy>=0
                        dWdx = (3*W(j,i)-4*W(j,i-1)+W(j,i-2))*0.5*hi; 
                        dV_dx = (3*V(j,i)-4*V(j,i-1)+V(j,i-2))*0.5*hi;
                    else
                        dWdx = -(3*W(j,i)-4*W(j,i+1)+W(j,i+2))*0.5*hi; 
                        dV_dx = -(3*V(j,i)-4*V(j,i+1)+V(j,i+2))*0.5*hi; 
                    end
                    if -dPsidx>=0
                        dWdy = (3*W(j,i)-4*W(j-1,i)+W(j-2,i))*0.5*hi; 
                        dV_dy = (3*V(j,i)-4*V(j-1,i)+V(j-2,i))*0.5*hi; 
                    else
                        dWdy = -(3*W(j,i)-4*W(j+1,i)+W(j+2,i))*0.5*hi;
                        dV_dy = -(3*V(j,i)-4*V(j+1,i)+V(j+2,i))*0.5*hi; 
                    end
                    LapW = (W(j-1,i)+W(j,i-1)-4*W(j,i)+W(j,i+1)+W(j+1,i))*hi^2;
                    LapV = (V(j-1,i)+V(j,i-1)-4*V(j,i)+V(j,i+1)+V(j+1,i))*hi^2;

                    % SF terms
                    LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;

                    % Determine new RHS
                    Sph(j,i) = (dPsidy*dWdx - dPsidx*dWdy + 2*V(j,i)*dVdy - eps*LapW) + Pres(j,i);
                    Swh(j,i) = (W(j,i) + LapPsi) + Wres(j,i);
                    Svh(j,i) = (dPsidy*dV_dx - dPsidx*dV_dy - eps*LapV) + Vres(j,i);
                else
                    % SF
                    Sph(j,Nx+1) = (Psi(j-1,Nx+1)-Psi(j-1,Nx)+Psi(j,Nx)-2*Psi(j,Nx+1)+Psi(j,Nx+2)+Psi(j+1,Nx+1)-Psi(j+1,Nx)) + Pres(j,Nx+1); 
                    
                    % V mom.
                    Svh(j,Nx+1) = (V(j-1,Nx+1)+V(j,Nx)-4*V(j,Nx+1)+V(j,Nx+2)+V(j+1,Nx+1)) + Vres(2,Nx+1);
                    
                end
            else
                if i<Nx+1
                    % V mom
                    if dPsidy>=0
                        dV_dx = (3*V(2,i)-4*V(2,i-1)+V(2,i-2))*0.5*hi; 
                    else
                        dV_dx = -(3*V(2,i)-4*V(2,i+1)+V(2,i+2))*0.5*hi; 
                    end
                    LapV = (V(1,i)+V(2,i-1)-4*V(2,i)+V(2,i+1)+V(3,i))*hi^2;
                    Svh(2,i) = (dPsidy*dV_dx - eps*LapV) + Vres(2,i);
                    
                else
                    Svh(2,Nx+1) = (V(1,Nx+1)+V(2,Nx)-4*V(2,Nx+1)+V(2,Nx+2)+V(3,Nx+1)) + Vres(2,Nx+1);
                end
            end
        end
    end
end