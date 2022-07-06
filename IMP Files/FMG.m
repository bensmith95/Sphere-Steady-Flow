%% FMG-FAS Alg. to solve the non-linear system
% Inputs: 
% Psi - Initial Stream Function (SF) guess     W - Initial Vorticity guess 
% V - Initial azimuthal velocity guess         
% BC - function that contains the Boundary Condition information
% H - size of fine grid                        Re - Reynolds Number
% Outputs:
% Psi - solution field for SF                  W - solution field for Vorticity
% V - solution field for azimuthal velocity    res - final residual of solution
function [Psi,W,V,res] = FMG(Psi,W,V,BC,H,Re)
    %% Pre-Condition and Restrict to Coarse Grid
    tol = 1e-5; % desired accuracy/tolerance 
    a = 0.9; % under-relaxation parameter
    h = H; str1 = fprintf('Restricting to Coarse Grid...\n');
    while h<2
        % Relax
        Sp = zeros(size(Psi)); Sw = Sp; Sv = Sp;
        [Psi,W,V] = Relax(Psi,W,V,Sp,Sw,Sv,BC,h,Re,9,a);
        
        % Restrict
        Psi = restrict(Psi); W = restrict(W); V = restrict(V);
        h = 2*h; tol = 4*tol;
        
        % Update BCs
        N = size(Psi); Nx = N(2)-2; Ny = N(1)-2; hi = 1/h; [Pxxt,Pt,VI] = BC(h,Re);
        % Psi
        Psi(Ny+1,2:Nx+1) = Pt(:)';
        % W
        W(Ny+1,2:Nx+1) = 2*(Psi(Ny+1,2:Nx+1)-Psi(Ny,2:Nx+1))*hi^2 - Pxxt(:)';
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(3:Ny,Nx+1) = (2*Psi(3:Ny,Nx+1)-Psi(2:Ny-1,Nx)-Psi(4:Ny+1,Nx))*hi^2;
        % V
        V(Ny+1,2:Nx+1) = VI;
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

    %% Solve on Coarse Grid
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving on %iX%i grid\n',[Nx,Ny]);
    Sp = zeros(size(Psi)); Sw = Sp; Sv = Sp;
    [Psi,W,V] = Relax(Psi,W,V,Sp,Sw,Sv,BC,h,Re,18,0.95);

    %% FMG Alg.
    while h>H
        % Interpolate 
        Psi = interp(Psi); W = interp(W); V = interp(V); h = h/2; tol = tol/4;
        % Update BCs
        N = size(Psi); Nx = N(2)-2; Ny = N(1)-2;
        fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving on %iX%i grid\n',[Nx,Ny]);
        hi = 1/h; [Pxxt,Pt,VI] = BC(h,Re);
        % Psi
        Psi(Ny+1,2:Nx+1) = Pt(:)';
        % W
        W(Ny+1,2:Nx+1) = 2*(Psi(Ny+1,2:Nx+1)-Psi(Ny,2:Nx+1))*hi^2 - Pxxt(:)';
        W(2:Ny+1,2) = -2*Psi(2:Ny+1,3)*hi^2;
        W(3:Ny,Nx+1) = (2*Psi(3:Ny,Nx+1)-Psi(2:Ny-1,Nx)-Psi(4:Ny+1,Nx))*hi^2;
        % V
        V(Ny+1,2:Nx+1) = VI;
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
        % Post-Relax
        Sp = zeros(size(Psi)); Sw = Sp; Sv = Sp;
        [Psi,W,V] = Relax(Psi,W,V,Sp,Sw,Sv,BC,h,Re,6,a);
        % FAS Alg.
        [Psi,W,V] = FAS(Psi,W,V,Sp,Sw,Sv,BC,h,Re,a,tol);
    end
    
    %% Determine Residual
    fprintf(repmat('\b',1,str1));
    [Pres,Wres,Vres] = residual(Psi,W,V,Sp,Sw,Sv,h,Re); 
    res = max([abs(Pres);abs(Wres);abs(Vres)],[],'all');
end