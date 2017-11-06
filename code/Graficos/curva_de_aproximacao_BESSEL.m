%% Calcula valores máximos de chirp

clear, clc, close all
F = linspace(1,10e9);
W = 2*pi*F;

D = 17e-6;
L = 100e3;
lambda = 1550e-9;
c = 299792458;

beta2 = -D*lambda^2/(2*pi*c);
beta2L = beta2*L;

% m = 16; % indíce de modulação global %
% Nc = 200;
mIM = 0.02;%sqrt(2*(m/100)^2./Nc); 
alfa = 5;
kappa = 15e12;
P0 = 5e-3;

phiIM = 0;
mFMc = alfa/2*(1 + kappa*P0./(1j*W))*1j*mIM*exp(1j*phiIM);
mFM = abs(mFMc);
phiFM = angle(mFMc);
dPSI = asin(1./sqrt(1 + (kappa*P0./W).^2));
% mFM = mIM*alfa/2*sqrt(1 + (kappa*P0./W).^2);
% dPSI = 
theta = -0.5*beta2L*2*pi*10e9*W;
u = 2*mFM.*sin(theta);

figure
plot(F, u)

% figure
% plot(F, besselj(u);
% xlabel('Frequência')
% ylabel('uk')
% 1;

% 
% figure
% plot(F, cos(theta), F, sin(theta))
% xlabel('Frequência')
% legend('cos(theta)', 'sin(theta)')
% 
% figure
% plot(theta, tan(theta), theta, theta)
% xlabel('theta')
% legend('tan')
% 
% figure
% plot(theta, cos(theta), theta, sin(theta))
% xlabel('theta')
% legend('cos(theta)', 'sin(theta)')

u = linspace(-0.6,0.6);
figure
plot(u, real(besselj(0,u)), '-b', 'LineWidth', 2);
hold on
plot(u, real(besselj(1,u)), '-k', 'LineWidth', 2);
plot(u, real(besselj(2,u)), '-g', 'LineWidth', 2);
plot(u, real(besselj(3,u)), '-r', 'LineWidth', 2);
plot(u, 1 - u.^2/4, '.b', 'MarkerSize', 10)
plot(u, u/2, '.k', 'MarkerSize', 10)
xlabel('u = 2m_{FM}sen(\theta)', 'FontSize', 14)
ylabel('J_n(u)', 'FontSize', 14)
legend('n = 0','n = 1', 'n = 2', 'n = 3', '1 - u^2/4', 'u/2')
set(gca, 'FontSize', 12)
axis([-0.6 0.6 -0.4 1])
set(gca, 'ytick', -0.4:0.2:1)
set(gca, 'xtick', -0.6:0.2:0.6)


% figure
% plot(u, real(besselj(0,u)), u, real(besselj(1,u)), u, real(besselj(2,u)), u, real(besselj(3,u)))
% xlabel('Frequência')
% ylabel('J(nk,uk)')
% legend('nk = 0','nk = 1', 'nk = 2', 'nk = 3')

% figure
% plot(F, kappa*P0./W)
% ylabel('adibático/transiente')

% figure
% plot(F, mFM)
% xlabel('Frequência')
% ylabel('mFM')

% figure, subplot(211)
% plot(F, rad2deg(dPSI))
% xlabel('Frequência')
% ylabel('deg')
% subplot(212)
% plot(dPSI, sin(dPSI), dPSI, cos(dPSI))
% legend('sin', 'cos')
% xlabel('Frequência')

fig = figure;
axes1 = axes('Parent',fig,...
    'Position',[0.13 0.11 0.599411764705882 0.815],...
    'Layer','top',...
    'FontSize',12);
box(axes1,'on');
hold(axes1,'all');
[f,kP] = meshgrid(linspace(1,10), linspace(1*10,5*15));
Z = 20*log10(kP./(2*pi*f));
[C,h] = contourf(f,kP,Z);
clabel(C,h);
xlabel('\Omega_k/2\pi (GHz)', 'FontSize', 14);
ylabel('\kappa P_0 (GHz)', 'FontSize', 14);
set(gca, 'FontSize', 12)
colorbar('peer',axes1,'FontSize',14);
set(axes1,'Position',[0.13 0.11 0.599411764705882 0.815]);
set(gca, 'xtick', 1:10)
% Create arrow
annotation(fig,'arrow',[0.843649312516885 0.843649312516885],...
    [0.735714285714286 0.911904761904762]);

% Create textbox
annotation(fig,'textbox',...
    [0.843537185956258 0.751219531880552 0.168115942028986 0.148305084745763],...
    'String',{'Chirp','Adiabático','Dominantnte'},...
    'FontSize',12,...
    'LineStyle','none');

% Create textbox
annotation(fig,'textbox',...
    [0.849175341258664 0.063125986429949 0.168115942028986 0.148305084745763],...
    'String',{'Chirp','Transiente','Dominantnte'},...
    'FontSize',12,...
    'LineStyle','none');

% Create arrow
annotation(fig,'arrow',[0.841781874039939 0.84331797235023],...
    [0.203389830508475 0.105932203389831]);


%% theta
fig = figure;
axes1 = axes('Parent',fig,...
    'Position',[0.13 0.11 0.599411764705882 0.815],...
    'Layer','top',...
    'FontSize',12);
[wimp,wk] = meshgrid(linspace(1,10), linspace(1,10));
Z = sin(-1/2*beta2*20e3*2*pi*2*pi*wimp.*wk*1e18);
% mesh(wk,wimp,Z)
[C,h] = contourf(wk,wimp,Z);
clabel(C,h);
xlabel('\Omega_k/2\pi (GHz)', 'FontSize', 14);
ylabel('\Omega_{IMP}/2\pi (GHz)', 'FontSize', 14);
set(gca, 'FontSize', 12)
colorbar('peer',axes1,'FontSize',14);
set(gca, 'xtick', 1:10)

fig = figure;
axes1 = axes('Parent',fig,...
    'Position',[0.13 0.11 0.599411764705882 0.815],...
    'Layer','top',...
    'FontSize',12);
[wimp,wk] = meshgrid(linspace(1,10), linspace(1,10));
Z = sin(-1/2*beta2*100e3*2*pi*2*pi*wimp.*wk*1e18);
% mesh(wk,wimp,Z)
[C,h] = contourf(wk,wimp,Z);
clabel(C,h);
xlabel('\Omega_k/2\pi (GHz)', 'FontSize', 14);
ylabel('\Omega_{IMP}/2\pi (GHz)', 'FontSize', 14);
set(gca, 'FontSize', 12)
colorbar('peer',axes1,'FontSize',14);
set(gca, 'xtick', 1:10)