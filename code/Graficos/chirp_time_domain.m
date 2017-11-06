%% Gera curvas de sinais com chirp no domínio do tempo

clear, clc, close all

t = linspace(-50e-9,50e-9,2^14);
T0 = 10e-9;
A = @(t,C,m) exp(-(1 + 1j*C)/2*(t/T0).^(2*m));
N = 2^14;

alfa = -5;
kappa =0.5;
P0 = 1e-3;

Pt = P0*(1 + A(t,0,2));
figure, plot(t, Pt, '-k', 'LineWidth', 4)
set(gca, 'xtick', []);
set(gca, 'ytick', []);
% xlabel('Tempo', 'FontSize', 14);
% ylabel('Potência Óptica', 'FontSize', 14);
set(gca, 'box', 'off')
set(gca, 'xcolor', [1 1 1])
set(gca, 'ycolor', [1 1 1])
axis([t(1) t(end) min(Pt)*0.9 1.2*max(Pt)])

dv = -alfa/(4*pi)*(diff([log(Pt(1)) log(Pt)]) + kappa*Pt);
figure, plot(t, dv, '-k', 'LineWidth', 4)
set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gca, 'box', 'off')
set(gca, 'xcolor', [1 1 1])
set(gca, 'ycolor', [1 1 1])
% xlabel('Tempo', 'FontSize', 14);
% ylabel('Desvio de Frequência', 'FontSize', 14);
axis([t(1) t(end) min(dv)*0.9 1.2*max(dv)]);

figure, plot(t, sin(2*pi/T0*t), '-k', 'LineWidth', 4);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gca, 'box', 'off')
set(gca, 'xcolor', [1 1 1])
set(gca, 'ycolor', [1 1 1])
% xlabel('Tempo', 'FontSize', 14);
% ylabel('Portadora Óptica', 'FontSize', 14);
% axis([t(1) t(end) min(dv)*0.9 1.2*max(dv)]);


figure, plot(t, sin(2*pi/T0*t + 12e3*dv), '-k', 'LineWidth', 4);
set(gca, 'xtick', []);
set(gca, 'ytick', []);
set(gca, 'box', 'off')
set(gca, 'xcolor', [1 1 1])
set(gca, 'ycolor', [1 1 1])
% xlabel('Tempo', 'FontSize', 14);
% ylabel('Portadora Óptica', 'FontSize', 14);
% axis([t(1) t(end) min(dv)*0.9 1.2*max(dv)]);