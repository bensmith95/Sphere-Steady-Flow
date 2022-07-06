%% Solves the radial jet boundary layer equations 
% Inputs:
% Re - square root of Reynolds number

function RadialJet(Re)
    %% Initialise
    % Load Inlet
    filename = ['../Flows/IMP/IMP_Re=',num2str(Re),'.mat'];
    load(filename,'VelIMP', 'beta', 'eta'); 
    UIMP = VelIMP{1}; VIMP = VelIMP{2}; WIMP = VelIMP{3}; Eta=eta;
    [~,ii] = min(abs(Eta-5)); % Determine etamax in eta 
    VBC = VIMP(end,ii:end); % V BC

    % Load BCs
    filename = '../Flows/BL.mat';
    load(filename,'VelBL','eta','theta'); VBL = VelBL{2};   
    Tm = -max(beta)/Re + pi/2; [~,j] = min(abs(theta-Tm)); % Determine betamax in theta 
    [~,i] = min(abs(eta-Eta(end))); % Determine Etamax in eta 
    VBC = [VBC VBL(j,i+1:end)]; % V BC

    % Space marching parameters
    Nbeta = length(beta); beta = -flip(beta);
    r = 1+[Eta(ii-1:end) eta(i+1:end)]/Re; dr = r(end)-r(end-1);
    r = [r logspace(log10(r(end)+dr),log10(20),200)]; Nr = length(r);  

    % Initialise flow field
    Urj = zeros(Nbeta,Nr); Vrj = Urj; Wrj = Urj; 

    % Inlet conditions
    Uin = -flip(UIMP(:,ii-1:ii))*sqrt(Re); Vin = flip(VIMP(:,ii-1:ii)); 
    Win = flip(WIMP(:,ii-1:ii))/sqrt(Re); Win(Win<0) = 0;
    Urj(:,1:2) = Uin; Wrj(:,1:2) = Win; 
    Vrj(:,1:2) = Vin; Vrj(1,2:length(VBC)+1) = VBC;  

    %% NEWTON-PICARD MULTIVARIATE SCHEME
    str1 = fprintf('Solving r =  %.3f\n',r(2));
    for i=3:Nr
        fprintf(repmat('\b',1,str1)); str1 = fprintf('Solving r =  %.3f\n',r(i));

        % Picard Method
        [Urj,Vrj,Wrj] = picard(Urj,Vrj,Wrj,r,beta,i,0.5);

        % Newton Method
        [Urj,Vrj,Wrj] = newton(Urj,Vrj,Wrj,r,beta,i,0.95);

    end
    Urj = Urj(:,2:end); Vrj = Vrj(:,2:end); Wrj = Wrj(:,2:end); r = r(2:end);
    fprintf(repmat('\b',1,str1)); 
    str1 = fprintf('All radial distances solved for.\n'); 

    % save file
    VelRJ{1} = Urj; VelRJ{2} = Vrj; VelRJ{3} = Wrj;
    filename = ['../Flows/RJ/RJ_Re=',num2str(Re),'.mat'];
    if exist('../Flows/RJ','dir')==0; mkdir('../Flows/RJ'); end
    save(filename,'VelRJ','r','beta');
    str2 = fprintf('Flow saved in %s\n', filename); pause(1)

    fprintf(repmat('\b',1,str2)); fprintf(repmat('\b',1,str1));
end