%% Solves the Impinging region equations 
% Inputs:
% Re - Reynolds Number (sqrt) 

function IMPRegion(Re)
    %% Initialise
    str1 = fprintf('Initialising...\n');
    % step size
    h = 2^-4;
    % set eta
    eta = 0:h:30; Neta = length(eta); 
    % set beta
    beta = 0:h:10; Nbeta = length(beta); 

    % assign variable storage
    Psi = zeros(Nbeta+2,Neta+2); 
    Omega = zeros(Nbeta+2,Neta+2); 
    V = zeros(Nbeta+2,Neta+2); 

    % Obtain and update BC's
    hi = 1/h; [Pxxt,Pt,VI] = BC(h,Re);
    % Psi
    Psi(Nbeta+1,2:Neta+1) = Pt(:)'; 
    % W
    Omega(Nbeta+1,2:Neta+1) = 2*(Psi(Nbeta+1,2:Neta+1)-Psi(Nbeta,2:Neta+1))*hi^2 - Pxxt(:)';
    Omega(2:Nbeta+1,2) = -2*Psi(2:Nbeta+1,3)*hi^2;
    Omega(3:Nbeta,Neta+1) = (2*Psi(3:Nbeta,Neta+1)-Psi(2:Nbeta-1,Neta)-Psi(4:Nbeta+1,Neta))*hi^2;
    % V
    V(Nbeta+1,2:Neta+1) = VI; V(2:Nbeta+1,2) = 1; 
    % Extrapolated points
    Psi(2:Nbeta+1,1) = Psi(2:Nbeta+1,3);
    Psi(2:Nbeta+1,Neta+2) = Psi(2:Nbeta+1,Neta);
    Psi(1,2:Neta+1) = 3*Psi(2,2:Neta+1) - 3*Psi(3,2:Neta+1) + Psi(4,2:Neta+1);
    Psi(Nbeta+2,2:Neta+1) = Psi(Nbeta,2:Neta+1);
    Omega(2:Nbeta+1,1) = 3*Omega(2:Nbeta+1,2) - 3*Omega(2:Nbeta+1,3) + Omega(2:Nbeta+1,4);
    Omega(2:Nbeta+1,Neta+2) = 3*Omega(2:Nbeta+1,Neta+1) - 3*Omega(2:Nbeta+1,Neta) + Omega(2:Nbeta+1,Neta-1);
    Omega(1,2:Neta+1) = -3*Omega(3,2:Neta+1) + Omega(4,2:Neta+1);
    Omega(Nbeta+2,2:Neta+1) = 3*Omega(Nbeta+1,2:Neta+1) - 3*Omega(Nbeta,2:Neta+1) + Omega(Nbeta-1,2:Neta+1);
    V(2:Nbeta+1,1) = 3*V(2:Nbeta+1,2) - 3*V(2:Nbeta+1,3) + V(2:Nbeta+1,4);
    V(2:Nbeta+1,Neta+2) = V(2:Nbeta+1,Neta);
    V(1,2:Neta+1) = V(3,2:Neta+1);
    V(Nbeta+2,2:Neta+1) = 3*V(Nbeta+1,2:Neta+1) - 3*V(Nbeta,2:Neta+1) + V(Nbeta-1,2:Neta+1);

    %% FMG-FAS Alg.
    fprintf(repmat('\b',1,str1)); str1 = fprintf('FMG-FAS Algorithm...\n');
    [Psi,Omega,V,~] = FMG(Psi,Omega,V,@BC,h,Re); 
    fprintf(repmat('\b',1,str1)); str1 = fprintf('Solution Converged.\n');

    %% POST PROCESS
    str2 = fprintf('Post-Processing...\n');
    % load Boundary Layer flow 
    filename = '../Flows/BL.mat';
    load(filename,'VelBL','eta','theta'); UB = VelBL{1}; Eta = 0:h:30; 
    % Determine U inlet
    Tm = -10/Re + pi/2; [~,i] = min(abs(theta-Tm));
    UI = sqrt(Re)*spline(eta,-UB(i,:),Eta); % U in
    % Determine W via finite differences of SF
    W = zeros(Nbeta,Neta); 
    W(2:Nbeta-1,:) = (Psi(4:Nbeta+1,2:Neta+1)-Psi(2:Nbeta-1,2:Neta+1))*0.5*hi;
    W(1,:) = -(3*Psi(2,2:Neta+1)-4*Psi(3,2:Neta+1)+Psi(4,2:Neta+1))*0.5*hi;
    W(Nbeta,:) = 0; W(:,1) = 0;  
    % Determine U via finite differences of SF
    U = zeros(Nbeta,Neta);
    U(:,2:Neta-1) = -(Psi(2:Nbeta+1,4:Neta+1)-Psi(2:Nbeta+1,2:Neta-1))*0.5*hi;
    U(:,Neta) = 0; U(:,1) = 0; U(Nbeta,:) = UI; U(1,:) = 0; 
    
    % Solve Poisson Equation for Pressure
    PB = VelBL{4}; PI = spline(eta,PB(i,:),Eta); % P in
    P = Pressure(U,V(2:end-1,2:end-1),W,PI,h); 

    % save data to file
    VelIMP{1} = U; VelIMP{2} = V(2:end-1,2:end-1); VelIMP{3} = W; 
    VelIMP{4} = Psi(2:end-1,2:end-1); VelIMP{5} = Omega(2:end-1,2:end-1); VelIMP{6} = P;
    filename = ['../Flows/IMP/IMP_Re=',num2str(Re),'.mat'];
    if exist('../Flows/IMP','dir')==0; mkdir ../Flows/IMP; end
    eta = Eta; save(filename, 'VelIMP', 'eta', 'beta')
    fprintf(repmat('\b',1,str2)); str2 = fprintf('Flow saved in %s\n', filename); pause(1)

    fprintf(repmat('\b',1,str2)); fprintf(repmat('\b',1,str1));
end
%% BOUNDARY CONDITIONS
function [Pxx,Psit,VI] = BC(h,Re)
    hi = 1/h;
    % load Boundary Layer flow 
    filename = '../Flows/BL.mat';
    load(filename,'VelBL','eta','theta'); 
    UB = VelBL{1}; VB = VelBL{2}; Eta = 0:h:30; 
    % Determine betamax in theta 
    Tm = -10/Re + pi/2; [~,i] = min(abs(theta-Tm));
    % Determine U inlet
    UI = sqrt(Re)*spline(eta,-UB(i,:),Eta);
    % Determine V inlet 
    VI = spline(eta,VB(i,:),Eta);
    % Vorticity BC's
    Pxx = zeros(size(Eta)); Pxx(2:end-1) = (UI(3:end)-UI(1:end-2))*0.5*hi;
    Pxx(1) = -(3*UI(1)-4*UI(2)+UI(3))*0.5*hi; Pxx(end) = (3*UI(end)-4*UI(end-1)+UI(end-2))*0.5*hi; Pxx=-Pxx;
    % SF BC's
    Psit = -cumtrapz(Eta,UI); % Psi Top
end