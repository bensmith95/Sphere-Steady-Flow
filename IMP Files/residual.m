%% Residual of the SF-Vort. equations 
% Inputs:
% Psi - Approximation for SF             W - Approximation for Vort.
% V - Approximation for azi. vel.        Sp - RHS of SF eqn
% Sw - RHS of Vort. eqn                  Sv - RHS of V eqn
% h - size of grid                       Re - Reynolds Number 
% Outputs:
% Pres - residual of SF eqn              
% Wres - residual of Vort. eqn 
% Vres - residual of V eqn               
function [Pres,Wres,Vres] = residual(Psi,W,V,Sp,Sw,Sv,h,Re)
    %% Initialise
    N = size(Psi); Nx = N(2)-2; Ny = N(1)-2;
    Pres = zeros(N); Wres = zeros(N); Vres = zeros(N);
    eps = 1/sqrt(Re); hi = 1/h;
    
    %% Iterate
    for j = Ny:-1:2
        for i = 3:Nx+1
            dPsidx = (Psi(j,i+1)-Psi(j,i-1))*0.5*hi; 
            dPsidy = (Psi(j+1,i)-Psi(j-1,i))*0.5*hi;
            dVdy = (V(j+1,i)-V(j-1,i))*0.5*hi; 
            if j>2
                if i<Nx+1
                    % SF & V mom. terms 
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
                    
                    % SF
                    Pres(j,i) = Sp(j,i) - (dPsidy*dWdx - dPsidx*dWdy + 2*V(j,i)*dVdy - eps*LapW);

                    % Vort.
                    LapPsi = (Psi(j-1,i)+Psi(j,i-1)-4*Psi(j,i)+Psi(j,i+1)+Psi(j+1,i))*hi^2;
                    Wres(j,i) = Sw(j,i) - (W(j,i) + LapPsi);

                    % V mom.
                    Vres(j,i) = Sv(j,i) - (dPsidy*dV_dx - dPsidx*dV_dy - eps*LapV);
                else
                    % SF
                    Pres(j,Nx+1) = Sp(j,Nx+1) - (Psi(j-1,Nx+1)-Psi(j-1,Nx)+Psi(j,Nx)-2*Psi(j,Nx+1)+Psi(j,Nx+2)+Psi(j+1,Nx+1)-Psi(j+1,Nx));
                    
                    % V mom.
                    Vres(j,Nx+1) = Sv(j,Nx+1) - (V(j-1,Nx+1)+V(j,Nx)-4*V(j,Nx+1)+V(j,Nx+2)+V(j+1,Nx+1)); 
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
                    Vres(2,i) = Sv(2,i) - (dPsidy*dV_dx - eps*LapV);     
                else
                    Vres(2,Nx+1) = Sv(2,Nx+1) - (V(1,Nx+1)+V(2,Nx)-4*V(2,Nx+1)+V(2,Nx+2)+V(3,Nx+1));
                end
            end
        end
    end
end