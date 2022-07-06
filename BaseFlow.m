%% Obtains and assembles the full flow 
% Inputs:
% Re - Renolds number (sqrt)
function BaseFlow(Re)
    fprintf('\nObtaining Base Flow:\n')
    %% Obtain Boundary Layer region 
    % Check if data aleady exists
    BLfile = 'Flows/BL.mat';
    if exist(BLfile,'file')
        str1 = fprintf('1) Boundary Layer Region File Found.\n'); pause(1)
    else
        % Solve 
        str1 = fprintf('1) Solving Boundary Layer Region...\n');
        cd 'BL Files'; BLRegion; cd ..
    end
    fprintf(repmat('\b',1,str1)); fprintf('1) Obtained Boundary Layer Region.\n')

    %% Obtain Impinging region
    % Check if data aleady exists
    IMPfile = ['Flows/IMP/IMP_Re=',num2str(Re),'.mat'];
    if exist(IMPfile,'file')
        str2 = fprintf('2) Impinging Region File Found.\n'); pause(1)
    else
        % Solve
        str2 = fprintf('2) Solving Impinging Region...\n');
        cd 'IMP Files'; IMPRegion(Re); cd ..
    end
    fprintf(repmat('\b',1,str2)); fprintf('2) Obtained Impinging Region.\n')

    %% Obtain Radial Jet region
    % Check if data aleady exists
    RJfile = ['Flows/RJ/RJ_Re=',num2str(Re),'.mat'];
    if exist(RJfile,'file')
        str3 = fprintf('3) Radial Jet Region File Found.\n'); pause(1)
    else
        % Solve
        str3 = fprintf('3) Solving Radial Jet Region...\n');
        cd 'RJ Files'; RadialJet(Re); cd ..
    end
    fprintf(repmat('\b',1,str3)); fprintf('3) Obtained Radial Jet Region.\n')
    
    %% Assemble Full Base Flow
    % Check if data aleady exists
    FULLfile = ['Flows/FULL/FULL_Re=',num2str(Re),'.mat'];
    if exist(FULLfile,'file')
        str4 = fprintf('4) Full Base Flow File Found.\n'); pause(1)
    else
        str4 = fprintf('4) Building Full Base Flow...\n');
        % load Boundary Layer region
        load(BLfile); UBL = VelBL{1}; VBL = VelBL{2}; WBL = VelBL{3}; BLeta = eta;
        % load Impinging region
        load(IMPfile); UIMP = -flip(VelIMP{1}); VIMP = flip(VelIMP{2}); WIMP = flip(VelIMP{3}); IMPeta = eta;
        % load Radial Jet region
        beta = []; load(RJfile); URJ = VelRJ{1}; VRJ = VelRJ{2}; WRJ = VelRJ{3};
        
        % Sigmoid function
        [~,re] = min(abs(r-(1+IMPeta(end)/Re))); [~,ii] = min(abs(IMPeta-5));
        sigmoid = 1./(1+exp(-(r(1:re)-(1+10/Re))*Re));
        for k=1:length(sigmoid)
            Wsig(:,k) = (1-sigmoid(k))*WIMP(:,ii+k-1)/sqrt(Re) + sigmoid(k)*WRJ(:,k);
            Usig(:,k) = (1-sigmoid(k))*UIMP(:,ii+k-1)/sqrt(Re) + sigmoid(k)*URJ(:,k)/Re;
            Vsig(:,k) = (1-sigmoid(k))*VIMP(:,ii+k-1) + sigmoid(k)*VRJ(:,k);
        end
        WRJ = [WIMP(:,1:ii-1)/sqrt(Re) Wsig WRJ(:,re+1:end)];
        URJ = [UIMP(:,1:ii-1)/sqrt(Re) Usig URJ(:,re+1:end)/Re];
        VRJ = [VIMP(:,1:ii-1) Vsig VRJ(:,re+1:end)];
        r = [IMPeta(1:ii)/Re+1 r(2:end)];

        % interpolate to same r
        [~,idx] = min(abs((r-1)*Re-30)); % determine etamax in r
        UBL = spline(BLeta,UBL,(r(1:idx)-1)*Re); 
        VBL = spline(BLeta,VBL,(r(1:idx)-1)*Re);
        WBL = spline(BLeta,WBL,(r(1:idx)-1)*Re); 
        % determine beta->Inf in theta
        thetaIN = -10/Re + pi/2; [~,j] = min(abs(theta-thetaIN)); 
        theta = [theta(1:j-1) beta/Re+pi/2];   
        % attach RJ region to BL region
        UFULL = zeros(length(theta),length(r)); VFULL = UFULL; WFULL = UFULL;
        UFULL(1:j-1,1:idx) = UBL(1:j-1,:); UFULL(j:end,:) = URJ; 
        VFULL(1:j-1,1:idx) = VBL(1:j-1,:); VFULL(j:end,:) = VRJ;
        WFULL(1:j-1,1:idx) = WBL(1:j-1,:)/Re; WFULL(j:end,:) = WRJ;
        % save data
        VelFULL{1} = UFULL; VelFULL{2} = VFULL; VelFULL{3} = WFULL;
        if exist('Flows/FULL','dir')==0; mkdir('Flows/FULL'); end
        save(FULLfile,'VelFULL','r','theta');
        fprintf(repmat('\b',1,str4)); str4 = fprintf('Flow saved in %s\n', FULLfile); pause(1) 
    end
    fprintf(repmat('\b',1,str4)); fprintf('4) Obtained Full Base Flow.\n\n')
end