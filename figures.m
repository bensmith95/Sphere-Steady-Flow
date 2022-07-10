%% PLOT FIGURES

clear all; clc; close all

% Reynolds number
Re = 100;

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
title(t1,'Boundary Layer Flow','interpreter','latex','fontsize',12);

% plot U
h1(1) = nexttile(t1); hold on; set(gca,'Color','none'); x = 0:1e-3:1; area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,UBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); 
set(get(cb,'label'),'string','$U$','interpreter','latex', 'fontsize',11);

% plot V
h1(2) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,VBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
set(get(cb,'label'),'string','$V$','interpreter','latex', 'fontsize',11);

% plot W
h1(3) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,WBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([-0.9 1.6]); 
set(get(cb,'label'),'string','$\overline{W}$','interpreter','latex', 'fontsize',11);

% plot P
h1(4) = nexttile(t1); hold on; set(gca,'Color','none'); area(x,sqrt(1-x.^2),'FaceColor','white')
fn = pcolor(XBL,YBL,PBL); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); 
set(get(cb,'label'),'string','$\bar{P}$','interpreter','latex', 'fontsize',11);

% settings
set(h1,'Colormap', jet)
t1.TileSpacing = 'tight'; t1.Padding = 'compact'; 
set(gcf, 'Position',  [325, 55, 825, 725])

% labelling
grey = [0.4 0.4 0.4]; ticks = [1 1+7.5/50 1+15/50 1+22.5/50 1+30/50]; 
for i=1:4
    nexttile(i); 
    % axis ticks
    set(gca,'Xtick',[]); 
    yticks(ticks); yticklabels({'0','7.5','15','22.5','30'}); 
    yh = ylabel('$\eta$','interpreter','latex', 'fontsize',11); 
    yh.Position(1)=-0.15; yh.Position(2)=1+13/50; yh.Rotation=0;
    % grid
    plot([0 0], [0 ticks(end)], 'color', grey);
    plot([0 ticks(end)], [0 0], 'color', grey);
    for j=ticks
        x = 0:1e-3:j; plot(x,sqrt(j^2-x.^2),'color', grey)
    end
    for k=15:15:75
        plot([cos(k*pi/180) ticks(end)*cos(k*pi/180)], [sin(k*pi/180) ticks(end)*sin(k*pi/180)], 'color', grey)
        text(0.96*sin(k*pi/180), 0.96*cos(k*pi/180), [num2str(k),'^o'],'fontsize',10,'HorizontalAlignment','right')
    end
end

% save figs
%figname = ['BLRegion'.png']; saveas(figure(1),figname);

%% Plot Impinging Region

% Load flow
IMPfile = ['Flows/IMP/IMP_Re=',num2str(Re),'.mat'];
load(IMPfile); UIMP = VelIMP{1}; VIMP = VelIMP{2}; WIMP = VelIMP{3}; 
Psi = VelIMP{4}; Omega = VelIMP{5}; PIMP = VelIMP{6}; IMPeta = eta;

% Plot U-W Velocity Magnitude & vector field
figure(2); TT = ['Velocity Field: $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(IMPeta,beta,sqrt(WIMP.^2+UIMP.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
fc = colorbar; colormap('jet');
set(get(fc,'label'),'string','$\sqrt{W^2+U^2}$','interpreter','latex');
hold on; ii = 32; jj = 16; 
quiver(IMPeta(1:ii:end),beta(1:jj:end),WIMP(1:jj:end,1:ii:end),UIMP(1:jj:end,1:ii:end),'color','k'); 
xlabel('\eta','fontsize',12); ylabel('-\beta','rotation',0,'fontsize',12); title(TT,'interpreter','latex','fontsize',12);
set(gcf, 'Position',  [200, 200, 1200, 400]); 
%figfile = ['Vel_Re=',num2str(Re),'.png']; saveas(figure(2),figfile); pause(1)

% Plot Stream Function
figure(3); TT = ['$\psi(\eta,\beta)$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re)];
contourf(IMPeta,beta,Psi,15,'LineColor','none' ); 
colorbar; colormap('jet');
xlabel('\eta','fontsize',12); ylabel('-\beta','rotation',0,'fontsize',12); title(TT,'interpreter','latex','fontsize',12);
axis equal; set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['SF_Re=',num2str(Re),'.png']; saveas(figure(3),figfile); pause(1)

% Plot Vorticity
figure(4); TT = ['$|\omega(\eta,\beta)|$ Contour Plot: $\sqrt{R_e} =$ ',num2str(Re)];
contourf(IMPeta,beta,abs(Omega),logspace(-8,0,40),'LineColor','none'); 
colorbar; colormap('jet');
xlabel('\eta','fontsize',12); ylabel('-\beta','fontsize',12,'rotation',0); title(TT,'interpreter','latex','fontsize',12);
axis equal; set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['Vort_Re=',num2str(Re),'.png']; saveas(figure(4),figfile); pause(1)

% Plot V velocity component 
figure(5); TT = ['$V(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(IMPeta,beta,VIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar; colormap('jet');
xlabel('\eta','fontsize',12); ylabel('-\beta','fontsize',12,'rotation',0); title(TT,'interpreter','latex','fontsize',12);
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['V_Re=',num2str(Re),'.png']; saveas(figure(5),figfile); pause(1)

% Plot W velocity component 
figure(6); TT = ['$W(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
Wp = WIMP; Wp(Wp<0)=NaN; % remove negative/reverse flow
fn = pcolor(IMPeta,beta,Wp); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar; colormap('jet'); 
xlabel('\eta','fontsize',12); ylabel('-\beta','fontsize',12,'rotation',0); title(TT,'interpreter','latex','fontsize',12);
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['W_Re=',num2str(Re),'.png']; saveas(figure(6),figfile); pause(1)

% Plot U velocity component 
figure(7); TT = ['$U(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(IMPeta,beta,-UIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar; colormap('jet'); 
xlabel('\eta','fontsize',12); ylabel('-\beta','fontsize',12,'rotation',0); title(TT,'interpreter','latex','fontsize',12);
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['U_Re=',num2str(Re),'.png']; saveas(figure(7),figfile); pause(1)

% Plot Pressure 
figure(8); TT = ['$\bar{P}(\eta,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)];
fn = pcolor(eta,beta,PIMP); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colorbar; colormap('jet'); 
xlabel('\eta','fontsize',12); xlim([0 30])
ylabel('-\beta','fontsize',12,'rotation',0); ylim([0 10]) 
title(TT,'interpreter','latex','fontsize',12);
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['P_Re=',num2str(Re),'.png']; saveas(figure(8),figfile); pause(1)

%% Plot radial jet figures

% load radial jet
RJfile = ['Flows/RJ/RJ_Re=',num2str(Re),'.mat'];
load(RJfile); Urj = VelRJ{1}; Vrj = VelRJ{2}; Wrj = VelRJ{3}; r1 = r;

% Plot W component
figure(9); TT = ['$W_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(r1,beta,Wrj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar;
xlabel('$r$','interpreter','latex','fontsize',12); xlim([r(1) 10]);
ylabel('\beta','fontsize',12,'rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex','fontsize',12)
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['RJ_W_Re=',num2str(Re),'.png']; saveas(figure(9),figfile); pause(1)

% Plot U component
figure(10); TT = ['$\overline{U}_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(r1,beta,Urj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar; caxis([0 2])
xlabel('$r$','interpreter','latex','fontsize',12); xlim([r(1) 5]);
ylabel('\beta','fontsize',12,'rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex','fontsize',12)
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['RJ_U_Re=',num2str(Re),'.png']; saveas(figure(10),figfile); pause(1)

% Plot V component
figure(11); TT = ['$V_{RJ}(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(r1,beta,Vrj); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar;
xlabel('$r$','interpreter','latex','fontsize',12); xlim([r(1) 2]);
ylabel('\beta','fontsize',12,'rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex','fontsize',12)
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['RJ_V_Re=',num2str(Re),'.png']; saveas(figure(11),figfile); pause(1)

%% plot r-> InF Solutions
figure(12); t = tiledlayout(1,3); 
TT = ['$r\rightarrow\infty$ solutions at $\sqrt{R_e}=$',num2str(Re)];

% Calculate relevant constants
M = trapz(beta,r1(end)^2*Wrj(:,end).^2); N = trapz(beta,r1(end)^3*Wrj(:,end).*Vrj(:,end));
delta = sqrt(2)/(M*3)^(1/3); epsilon = N/M; 

% Plot U(r->inf,beta)
nexttile(t); plot(beta,Urj(:,end)); hold on;
Uinf = -sqrt(2)/(delta*r1(end))*tanh(beta/(sqrt(2)*delta));
plot(beta(1:4:end),Uinf(1:4:end),'kx')
legend('Numerical','Analytical','Location','southwest'); 
title('(a)','interpreter','latex','fontsize',10)
xlabel('\beta','fontsize',12); 
ylabel('$\overline{U}_\infty(\beta)$','interpreter','latex','fontsize',12);

% Plot V(r->inf,beta)
nexttile(t); plot(beta,Vrj(:,end)); hold on;
Vinf = epsilon/(delta^2*r1(end)^2)*sech(beta/(sqrt(2)*delta)).^2;
plot(beta(1:4:end),Vinf(1:4:end),'kx')
legend('Numerical','Analytical','Location','northwest')
title('(b)','interpreter','latex','fontsize',10)
xlabel('\beta','fontsize',12);
ylabel('$V_\infty(\beta)$','interpreter','latex','fontsize',12);

% Plot W(r->inf,beta)
nexttile(t); plot(beta,Wrj(:,end)); hold on;
Winf = 1/(delta^2*r1(end))*sech(beta/(sqrt(2)*delta)).^2;
plot(beta(1:4:end),Winf(1:4:end),'kx')
legend('Numerical','Analytical','Location','northwest')
title('(c)','interpreter','latex','fontsize',10)
xlabel('\beta','fontsize',12); 
ylabel('$W_\infty(\beta)$','interpreter','latex','fontsize',12);

set(gcf, 'Position',  [50, 200, 1450, 400]); t.TileSpacing = 'tight'; t.Padding = 'tight'; 
title(t,TT,'interpreter','latex','fontsize',12);
%figfile = ['RJinf_Re=',num2str(Re),'.png']; saveas(figure(12),figfile); pause(1)

%% Combine IMPINGING REGION & RADIAL JET
beta = -flip(beta);
UIMP = -flip(VelIMP{1}); VIMP = flip(VelIMP{2}); WIMP = flip(VelIMP{3});

% attach IMP region to RJ region
%[~,ii] = min(abs(eta-5)); R = [eta(1:ii)/Re+1 r(2:end)]; 
%WRJ = [WIMP(:,1:ii)/sqrt(Re) Wrj(:,2:end)]; 
%URJ = [UIMP(:,1:ii)/sqrt(Re) Urj(:,2:end)/Re];
%VRJ = [VIMP(:,1:ii) Vrj(:,2:end)];

% Define Sigmoid function
[~,re] = min(abs(r1-(1+IMPeta(end)/Re))); [~,ii] = min(abs(IMPeta-5));
sigmoid = 1./(1+exp(-(r1(1:re)-(1+10/Re))*Re));

% Combine solutions
for k=1:length(sigmoid)
    Wsig(:,k) = (1-sigmoid(k))*WIMP(:,ii+k-1)/sqrt(Re) + sigmoid(k)*Wrj(:,k);
    Usig(:,k) = (1-sigmoid(k))*UIMP(:,ii+k-1)/sqrt(Re) + sigmoid(k)*Urj(:,k)/Re;
    Vsig(:,k) = (1-sigmoid(k))*VIMP(:,ii+k-1) + sigmoid(k)*Vrj(:,k);
end
WRJ = [WIMP(:,1:ii-1)/sqrt(Re) Wsig Wrj(:,re+1:end)];
URJ = [UIMP(:,1:ii-1)/sqrt(Re) Usig Urj(:,re+1:end)/Re];
VRJ = [VIMP(:,1:ii-1) Vsig Vrj(:,re+1:end)];
R = [IMPeta(1:ii)/Re+1 r1(2:end)];

% Plot W component
figure(13); TT = ['$W(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(R,beta,WRJ); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar;
xlabel('$r$','interpreter','latex','fontsize',12); xlim([1 1.4]);
ylabel('\beta','fontsize',12,'rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex','fontsize',12)
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['WRJ_Re=',num2str(Re),'.png']; saveas(figure(13),figfile); pause(1)

% Plot U component
figure(14); TT = ['$U(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(R,beta,URJ); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar;
xlabel('$r$','interpreter','latex','fontsize',12); xlim([1 1.4]);
ylabel('\beta','fontsize',12,'rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex','fontsize',12)
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['URJ_Re=',num2str(Re),'.png']; saveas(figure(14),figfile); pause(1)

% Plot V component
figure(15); TT = ['$V(r,\beta)$ at $\sqrt{R_e} =$ ',num2str(Re)]; 
fn = pcolor(R,beta,VRJ); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
colormap('jet'); colorbar;
xlabel('$r$','interpreter','latex','fontsize',12); xlim([1 1.4]);
ylabel('\beta','fontsize',12,'rotation',0); set(gca, 'YDir','reverse')
title(TT,'interpreter','latex','fontsize',12)
set(gcf, 'Position',  [200, 200, 1200, 400])
%figfile = ['VRJ_Re=',num2str(Re),'.png']; saveas(figure(15),figfile); pause(1)

%% Plot full flow figures

% load full flow
FULLfile = ['Flows/FULL/FULL_Re=',num2str(Re),'.mat'];
load(FULLfile); U = VelFULL{1}; V = VelFULL{2}; W = VelFULL{3};

% transform to cartesian coords
X = repmat(sin(theta)',1,length(r)); Y = repmat(cos(theta)',1,length(r));
for i=1:length(r); X(:,i) = r(i)*X(:,i); Y(:,i) = r(i)*Y(:,i); end

% Plot velocity components
figure(16); t1 = tiledlayout(1,3); 
title(t1,['Velocity Components at $\sqrt{R_e}$=',num2str(Re)],'interpreter','latex','fontsize',12);

% plot U
h(1) = nexttile(t1);
fn = pcolor(X,Y,U); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); 
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','$U$','interpreter','latex', 'fontsize',11);

% plot V
h(2) = nexttile(t1); 
fn = pcolor(X,Y,V); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','$V$','interpreter','latex', 'fontsize',11);

% plot W
h(3) = nexttile(t1); 
fn = pcolor(X,Y,W); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 0.3]) 
xlim([0 1.25]); ylim([0 1.25]);
set(get(cb,'label'),'string','$W$','interpreter','latex', 'fontsize',11);

% settings
set(h,'Colormap', jet)
t1.TileSpacing = 'compact'; t1.Padding = 'compact'; 
set(gcf, 'Position',  [50, 150, 1450, 500])

%figname = ['FullComps_Re=',num2str(Re),'.png']; saveas(figure(16),figname);

%% Plot Velocity magnitudes and vector field
figure(17); t2 = tiledlayout(1,2); 
title(t2,['Flow at $\sqrt{R_e}=$',num2str(Re)],'interpreter','latex', 'fontsize',12)

% plot planar magnitude
h(1) = nexttile(t2);
fn = pcolor(X,Y,sqrt(U.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none');
cb = colorbar('southoutside'); 
xlim([0 2]); ylim([0 1.25]);
set(get(cb,'label'),'string','$\sqrt{U^2+W^2}$','interpreter','latex', 'fontsize',11);

% plot velcoity magnitude and vector field
h(2) = nexttile(t2); 
fn = pcolor(X,Y,sqrt(U.^2+V.^2+W.^2)); fn.FaceColor = 'interp'; set(fn,'EdgeColor','none'); 
cb = colorbar('southoutside'); caxis([0 1]); set(cb, 'YTick', 0:0.2:1);
xlim([0 1.2]); ylim([0 1.05]);
set(get(cb,'label'),'string','$||${\boldmath$u$}$||$','interpreter','latex', 'fontsize',11);
hold on;
% convert to cartesian coords
u = W.*sin(theta') + U.*cos(theta'); v = W.*cos(theta') - U.*sin(theta');
quiver(X(1:40:end,1:40:end),Y(1:40:end,1:40:end),u(1:40:end,1:40:end),v(1:40:end,1:40:end),'color',[0.7 0.7 0.7],'autoscalefactor',0.1)

% settings
set(h,'Colormap', jet)
t2.TileSpacing = 'compact'; t2.Padding = 'compact'; 
set(gcf, 'Position',  [50, 150, 1450, 500])

%figname = ['FullMags_Re=',num2str(Re),'.png']; saveas(figure(17),figname);