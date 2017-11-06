clear, close all

warning off

addpath algoritmos\ -end

%% Funções auxiliares
J = @(n,x) besselj(n,x);
propaga =@(Ein,L,Dhat) ifft(exp(L*Dhat).*fft(Ein));
binomialc = @(n,k) nchoosek(n,k);
t2f = @(x) fftshift(fft(x))/length(x);
corta = @(x, indices) x(indices);
P = @(x)  sum(abs(x).^2)/length(x);
erro = @(xteorico, x) 100*P(x - xteorico)/P(xteorico);
dc2ac = @(x) x - mean(x);

%% Parametros de entrada
%% OFDM
Nc = 128;           % Número de subportadoras
fc = 5e9;         % Frequência central (Hz)
BW = Nc*0.05e9;           % Bandwidth (Hz)
Deltaf = BW/Nc;     

if mod(Nc,2) == 0
    F = (-Nc/2:Nc/2-1)*Deltaf + fc; % Nc par
else
    F =  (-floor(Nc/2):floor(Nc/2))*Deltaf + fc;
end

W = 2*pi*F;

%% laser
mIM = 0.02*ones(size(F));
QPSK = [pi/4 3*pi/4 -3*pi/4 -pi/4];
PSK8 = [0 pi/4 pi/2 3*pi/4 pi -3*pi/4 -pi/2 -pi/4];
phiIM = QPSK(randi(3, [1, Nc]) + 1);

lambda = 1550e-9;
c = 299792458;

% Artigo Tang (30Gb/s Signal transmission over 40km DML...)
laserfreq = c/lambda;
alfa = 1;   
OpticalConfinementEff = 0.07;
GainCompressionFac = 7.5e-23;
Vact = 29.2e-18;
DifQuantumEffpFacet = 0.139;
hPlanck = 6.626068e-34;
% kappa = 2*OpticalConfinementEff*GainCompressionFac/(Vact*DifQuantumEffpFacet*hPlanck*laserfreq);
kappa = 1e13;  % 	 
P0 = 4e-3;

%% Dados para FFT
N2 = 1;
deltaf = Deltaf/N2;
N = 2^14;

if floor(F(1)/deltaf) ~= ceil(F(1)/deltaf)
    disp('Amostragem da FFT não está alinhada com espectro OFDM');
    pause
end    

fa = N*deltaf;
deltat = 1/fa;
t = 0:deltat:(N-1)*deltat;
f = -fa/2:deltaf:fa/2-deltaf;
indices = 1:N2:length(f);

%% fibra
D = 17e-6;
L = 40e3;

beta2 = -D*lambda^2/(2*pi*c);
theta = -0.5*beta2*W.^2*L;

%%
Pt = P0 + zeros(size(t));
x = zeros(size(t));
At = sqrt(P0);

mFM = zeros(size(W));
phiFM = zeros(size(W));
dPSI = zeros(size(W));

for k = 1:Nc
    FM = alfa/2*(1 + kappa*P0./(1j*W(k)))*1j*mIM(k)*exp(1j*phiIM(k));
    mFM(k) = abs(FM);
    phiFM(k) = angle(FM);
    dPSI(k) = phiFM(k) - phiIM(k);
    Pt = Pt + P0*mIM(k)*cos(W(k)*t+phiIM(k));
    x = x + kappa*P0*mIM(k)/W(k)*sin(W(k)*t+phiIM(k));
    At = At + sqrt(P0)*mIM(k)/2*cos(W(k)*t + phiIM(k));
end

if min(Pt) <= 0
    disp('erro! Potência negativa!!!!!!!!!!!!!!!!!!')
    pause
end

fprintf('Indíce de Modulação Geral(+): %.3f\n', max(Pt)/mean(Pt) - 1);
fprintf('Indíce de Modulação Geral(-): %.3f\n', 1 - min(Pt)/mean(Pt));
fprintf('ER = %f [dB]\n', 10*log10(max(Pt)/min(Pt)))

%% Calcula Modelos
%% Modelo teórico
disp('Modelo Teórico')
TeoricoNumerico
teorico.I = dc2ac(teorico.I);
teorico.Idet = t2f(teorico.I);

%% 
figure
f = f/1e9;
plot(corta(f, indices), corta(abs(teorico.Idet), indices), 'ok')
hold on
plot(corta(f, indices), corta(abs(t2f(abs(teorico.E0).^2)),indices), 'sr')
legend(sprintf('Teórico Z = %.1f km', L/1e3), 'Teórico Z = 0 km')
axis([0 2*max(F)/1e9 0 0.1e-3])
