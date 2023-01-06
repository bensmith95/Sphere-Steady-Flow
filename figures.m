%% PLOT FIGURES

clear all; clc; close all

% Reynolds number
Re = sqrt(1e4);

% display info
set(0,'units','pixels'); disp = get(0,'ScreenSize');
Lx = disp(3); Ly = disp(4); 

%% Plot Boundary Layer Region
% Note - not to scale

% load file
BLfile = 'Flows/BL.mat'; load(BLfile); 
UBL = VelBL{1}; VBL = VelBL{2}; WBL = VelBL{3}; PBL = VelBL{4};
BLeta = eta; BLtheta = theta;

% initialise co-ords
rBL = 1+BLeta/50;
XBL = repmat(sin(BLtheta)',1,length(rBL)); YBL = repmat(cos(BLtheta)',1,length(rBL));
for i=1:length(rBL); XBL(:,i) = rBL(i)*XBL(:,i); YBL(:,i) = rBL(i)*YBL(:,i); end

figure(1); t1 = tiledlayout(2,2); 
title(t1,'Boundary Layer Flow','interpreter','latex');

% plot U
h1(1) = nexttile(t1); hold on; set(gca,'Color','none'); x = 0:1e-3:1; area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,UBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); caxis([0 0.15]);
set(get(cb,'label'),'string','(a) $U$','interpreter','latex');

% plot V
h1(2) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,VBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
set(get(cb,'label'),'string','(b) $V$','interpreter','latex');

% plot W
h1(3) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,WBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([-1 2.2]); 
set(get(cb,'label'),'string','(c) $\overline{W}$','interpreter','latex');

% plot P
h1(4) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,PBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([-3 0]);
set(get(cb,'label'),'string','(d) $\bar{P}$','interpreter','latex');

% settings
set(h1,'Colormap', jet)
t1.TileSpacing = 'tight'; t1.Padding = 'compact'; 
set(gcf, 'Position',  [0.225*Lx, 0.075*Ly, 0.55*Lx, 0.825*Ly])

% labelling
grey = [0.4 0.4 0.4]; ticks = [1 1+7.5/50 1+15/50 1+22.5/50 1+30/50]; 
for i=1:4
    nexttile(i); 
    % axis ticks
    set(gca,'Xtick',[]); 
    yticks(ticks); yticklabels({'0','7.5','15','22.5','30'}); 
    yh = ylabel('$\eta$','interpreter','latex'); 
    yh.Position(1)=-0.15; yh.Position(2)=1+13/50; yh.Rotation=0;
    % grid
    plot([0 0], [0 ticks(end)], 'color', grey);
    plot([0 ticks(end)], [0 0], 'color', grey);
    for j=ticks
        x = 0:1e-3:j; plot(x,sqrt(j^2-x.^2),'color', grey)
    end
    for k=15:15:75
        plot([cos(k*pi/180) ticks(end)*cos(k*pi/180)], [sin(k*pi/180) ticks(end)*sin(k*pi/180)], 'color', grey)
        text(0.96*sin(k*pi/180), 0.96*cos(k*pi/180), [num2str(k),'^o'],'HorizontalAlignment','right')
    end
end

%% Plot Impinging Region

% Load flow
IMPfile = ['Flows/IMP/IMP_Re=',num2str(Re),'.mat'];
load(IMPfile); UIMP = VelIMP{1}; VIMP = VelIMP{2}; WIMP = VelIMP{3}; 
Psi = VelIMP{4}; Omega = VelIMP{5}; PIMP = VelIMP{6}; IMPeta = eta; IMPbeta = beta;

% Plot U-W Velocity Magnitude & vector field
figure(2); TT = ['Impinging Region Velocity Field: $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(IMPeta,beta,sqrt(WIMP.^2+UIMP.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
fc = colorbar('location','southoutside'); caxis([0,0.35]); colormap('jet');
set(get(fc,'label'),'string','$\sqrt{W^2+U^2}$','interpreter','latex');
hold on; ii = 16; jj = 16; 
quiver(IMPeta(1:ii:end),beta(1:jj:end),WIMP(1:jj:end,1:ii:end),UIMP(1:jj:end,1:ii:end),'color','k'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly]); 

% Plot Stream Function
figure(3); TT = ['$\psi(\eta,\beta)$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re)];
contourf(IMPeta,beta,Psi,15,'LineColor','none' ); 
colorbar('location','southoutside'); caxis([0,0.8]); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

% Plot Vorticity
figure(4); TT = ['$\omega(\eta,\beta)$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re)];
contourf(IMPeta,beta,Omega,15,'LineColor','none'); 
colorbar('location','southoutside'); caxis([-0.15,0.1]); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

% Plot V velocity component 
figure(5); TT = ['$V_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(IMPeta,beta,VIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); colormap('jet');
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

% Plot W velocity component 
figure(6); TT = ['$W_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
Wp = WIMP; Wp(Wp<0)=NaN; % remove negative/reverse flow
fn = pcolor(IMPeta,beta,Wp); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([0,0.35]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

% Plot U velocity component 
figure(7); TT = ['$U_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(IMPeta,beta,-UIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([0,0.15]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

% Plot Pressure 
figure(8); TT = ['$\bar{P}_{IMP}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(eta,beta,PIMP/sqrt(Re)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar('location','southoutside'); caxis([-0.015,0.01]); colormap('jet'); 
xlabel('\eta'); ylabel('-\beta','rotation',0); title(TT,'interpreter','latex');
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

%% Plot radial jet figures

% load radial jet
RJfile = ['Flows/RJ/RJ_Re=',num2str(Re),'.mat'];
load(RJfile); Urj = VelRJ{1}; Vrj = VelRJ{2}; Wrj = VelRJ{3}; r1 = r;

% Plot W component
figure(9); TT = ['$W_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(r1,beta,Wrj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar('location','southoutside'); caxis([0,0.35]);
xlabel('$r$','interpreter','latex'); xlim([r(1) 10]);
ylabel('\beta','rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex')
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

% Plot U component
figure(10); TT = ['$\overline{U}_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(r1,beta,Urj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar('location','southoutside'); caxis([0,3]);
xlabel('$r$','interpreter','latex'); xlim([r(1) 2]);
ylabel('\beta','rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex')
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

% Plot V component
figure(11); TT = ['$V_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(r1,beta,Vrj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar('location','southoutside'); caxis([0,1]);
xlabel('$r$','interpreter','latex'); xlim([r(1) 2]);
ylabel('\beta','rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex')
set(gcf, 'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

%% plot r -> InF Solutions
figure(12); t = tiledlayout(1,3); 
TT = ['$r\rightarrow\infty$ Radial Jet solutions at $\sqrt{R_e}=$ ',num2str(Re)];

% Calculate relevant constants
M = trapz(beta,r1(end)^2*Wrj(:,end).^2); N = trapz(beta,r1(end)^3*Wrj(:,end).*Vrj(:,end));
K = sqrt(2)/(M*3)^(1/3); sigma = N/M; 

% Plot U(r->inf,beta)
nexttile(t); plot(beta,Urj(:,end),'-b'); hold on;
Uinf = -sqrt(2)/(K*r1(end))*tanh(beta/(sqrt(2)*K));
plot(beta(1:4:end),Uinf(1:4:end),'kx')
legend('Numerical','Analytical','Location','southwest'); 
title('(a) $\overline{U}_\infty$','interpreter','latex')
xlabel('\beta'); 

% Plot V(r->inf,beta)
nexttile(t); plot(beta,Vrj(:,end),'-b'); hold on;
Vinf = sigma/(K^2*r1(end)^2)*sech(beta/(sqrt(2)*K)).^2;
plot(beta(1:4:end),Vinf(1:4:end),'kx')
legend('Numerical','Analytical','Location','northwest')
title('(b) $V_\infty$','interpreter','latex')
xlabel('\beta');

% Plot W(r->inf,beta)
nexttile(t); plot(beta,Wrj(:,end),'-b'); hold on;
Winf = 1/(K^2*r1(end))*sech(beta/(sqrt(2)*K)).^2;
plot(beta(1:4:end),Winf(1:4:end),'kx')
legend('Numerical','Analytical','Location','northwest')
title('(c) $W_\infty$','interpreter','latex') 
xlabel('\beta');

t.TileSpacing = 'compact'; t.Padding = 'tight';
title(t,TT,'interpreter','latex');
set(gcf,'Position',  [0.1*Lx, 0.2*Ly, 0.8*Lx, 0.6*Ly])

%% Plot full flow figures

% load full flow
FULLfile = ['Flows/FULL/FULL_Re=',num2str(Re),'.mat'];
load(FULLfile); U = VelFULL{1}; V = VelFULL{2}; W = VelFULL{3};

% transform to cartesian coords
X = repmat(sin(theta)',1,length(r)); Y = repmat(cos(theta)',1,length(r));
for i=1:length(r); X(:,i) = r(i)*X(:,i); Y(:,i) = r(i)*Y(:,i); end

% Plot velocity components
figure(13); t1 = tiledlayout(1,3); 
title(t1,['Velocity Components at $\sqrt{R_e}$= ',num2str(Re)],'interpreter','latex');

% plot U
h(1) = nexttile(t1);
fn = pcolor(X,Y,U); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); caxis([0 0.15]);
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','(a) $U$','interpreter','latex');

% plot V
h(2) = nexttile(t1); 
fn = pcolor(X,Y,V); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','(b) $V$','interpreter','latex');

% plot W
h(3) = nexttile(t1); 
fn = pcolor(X,Y,W); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 0.35]) 
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','(c) $W$','interpreter','latex');

% settings
set(h,'Colormap', jet)
t1.TileSpacing = 'compact'; t1.Padding = 'compact';
set(gcf,'Position',  [0.1*Lx, 0.25*Ly, 0.8*Lx, 0.5*Ly])

%% Plot Velocity magnitudes and vector field
figure(14); t2 = tiledlayout(1,2); 
title(t2,['Flow at $\sqrt{R_e}=$ ',num2str(Re)],'interpreter','latex')

% plot planar magnitude
h(1) = nexttile(t2);
fn = pcolor(X,Y,sqrt(U.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); 
xlim([0 2]); ylim([0 1.25]);
set(get(cb,'label'),'string','$\sqrt{U^2+W^2}$','interpreter','latex');

% plot velcoity magnitude and vector field
h(2) = nexttile(t2); 
fn = pcolor(X,Y,sqrt(U.^2+V.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([0 1.2]); ylim([0 1.05]);
set(get(cb,'label'),'string','$||${\boldmath$U$}$||$','interpreter','latex');
hold on;
% convert to cartesian coords
u = W.*sin(theta') + U.*cos(theta'); v = W.*cos(theta') - U.*sin(theta');
quiver(X(1:40:end,1:40:end),Y(1:40:end,1:40:end),u(1:40:end,1:40:end),v(1:40:end,1:40:end),'color',[0.7 0.7 0.7],'autoscalefactor',0.1)

% settings
colormap('jet') 
set(h,'Colormap', jet) 
t2.TileSpacing = 'compact'; t2.Padding = 'compact'; 
set(gcf,'Position',  [0.1*Lx, 0.25*Ly, 0.8*Lx, 0.5*Ly])

%% Plot Velocities near Equator
figure(15); t3 = tiledlayout(1,3); 
title(t3,['Equatorial Velocity Components at $\sqrt{R_e}$= ',num2str(Re)],'interpreter','latex');

% plot U
h(1) = nexttile(t3);
fn = pcolor(X,Y,U); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); caxis([0 0.15]);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','(a) $U$','interpreter','latex');

% plot V
h(2) = nexttile(t3); 
fn = pcolor(X,Y,V); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','(b) $V$','interpreter','latex');

% plot W
h(3) = nexttile(t3); 
fn = pcolor(X,Y,W); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 0.35]) 
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','(c) $W$','interpreter','latex');

% settings
set(h,'Colormap', jet)
t3.TileSpacing = 'compact'; t3.Padding = 'compact'; 
set(gcf, 'Position',  [0.1*Lx, 0.25*Ly, 0.8*Lx, 0.5*Ly])

%% Plot Velocity magnitudes and vector field near Equator 

figure(16); t4 = tiledlayout(1,2); 
title(t4,['Equatorial Flow at $\sqrt{R_e}=$',num2str(Re)],'interpreter','latex')

% plot planar magnitude
h(1) = nexttile(t4);
fn = pcolor(X,Y,sqrt(U.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); caxis([0 0.35]);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','$\sqrt{U^2+W^2}$','interpreter','latex');

% plot velcoity magnitude and vector field
h(2) = nexttile(t4); 
fn = pcolor(X,Y,sqrt(U.^2+V.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([1-10/Re 1+30/Re]); ylim([0 (IMPbeta(end)+2.5)/Re]);
set(get(cb,'label'),'string','$||${\boldmath$U$}$||$','interpreter','latex');
hold on;
quiver(X(1:20:end,1:20:end),Y(1:20:end,1:20:end),u(1:20:end,1:20:end),v(1:20:end,1:20:end),'color','black','autoscalefactor',0.05)

% settings
set(h,'Colormap', jet)
t4.TileSpacing = 'compact'; t4.Padding = 'compact'; 
set(gcf, 'Position',  [0.1*Lx, 0.25*Ly, 0.8*Lx, 0.5*Ly])
