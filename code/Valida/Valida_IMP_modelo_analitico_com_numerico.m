 clear, close all

addpath('M');

%% Funções auxiliares
J = @(n,x) besselj(n,x);
propaga =@(Ein,L,Dhat) ifft(exp(L*Dhat).*fft(Ein));
binomialc = @(n,k) nchoosek(n,k);
t2f = @(x) fftshift(fft(x))/length(x);
f2t = @(x) ifft(ifftshift(x))*length(x);
corta = @(x, indices) x(indices);
P = @(x)  sum(abs(x).^2)/length(x);
erro = @(xteorico, x) 100*P(x - xteorico)/P(xteorico);
dc2ac = @(x) x - mean(x);

%% Parametros de entrada
%% OFDM
Nc = 7;             % Número de subportadoras
f0 = 4e9;           % Frequência central (Hz)
Deltaf = 0.5e9;
BW = Nc*Deltaf;     % Bandwidth (Hz)

F = (0:Nc-1)*Deltaf + f0; % Nc par

if Nc*Deltaf > f0
    disp('Erro na faixa de frequência')
    pause
elseif F(1) <= 0 || F(end) > 10e9
    disp('Fora da faixa de frequências');
    pause
end

W = 2*pi*F;

%% laser
mIM = 0.05*ones(size(F));
QPSK = [pi/4 3*pi/4 -3*pi/4 -pi/4];
PSK8 = [0 pi/4 pi/2 3*pi/4 pi -3*pi/4 -pi/2 -pi/4];
phiIM = 0*QPSK(randi(4, [1, Nc]));
% phiIM = pi*ones(size(mIM));

lambda = 1550e-9;
c = 299792458;

alfa = 5;
kappa = 0*3e13;  % 	 
P0 = 4e-3;

%% Dados para FFT
N2 = 1;
deltaf = Deltaf/N2;
N = 2^14;
deltat = 1/(N*deltaf);
fa = N*deltaf;
t = 0:deltat:(N-1)*deltat;
f = -fa/2:deltaf:fa/2-deltaf;
indices = 1:N2:length(f);

if fa <= 2*F(end)
    disp('Erro Frequ: fa <= 2*F(end)');
    pause
end

for k = 1:length(F)
    if isempty(find(abs(f) == F(k), 1))
        disp('Amostragem da FFT não está alinhada com espectro OFDM');
        pause
    end
end    

%% fibra
D = 17e-6;
L = 100e3;

beta2 = -D*lambda^2/(2*pi*c);
beta2L = beta2*L;
Dhat = -1j*0.5*beta2*(2*pi*ifftshift(f)).^2;

%%
IM = 0;
AM = 0;
FM = 0;
Adiab = 0;

mFM = zeros(size(W));
phiFM = zeros(size(W));
dPSI = zeros(size(W));

for k = 1:Nc
    mFMc = alfa/2*(1 + kappa*P0./(1j*W(k)))*1j*mIM(k)*exp(1j*phiIM(k));
    mFM(k) = abs(mFMc);
    phiFM(k) = angle(mFMc);
    dPSI(k) = phiFM(k) - phiIM(k);
    Adiab = Adiab + kappa*P0*mIM(k)/W(k)*sin(W(k)*t+phiIM(k));
    AM = AM + mIM(k)/2*cos(W(k)*t + phiIM(k));
    IM = IM + mIM(k)*cos(W(k)*t + phiIM(k));
    FM = FM + mFM(k)*sin(W(k)*t + phiFM(k));   
end

Pt = P0*(1 + IM);

if min(Pt) <= 0
    disp('Erro! Potência negativa!!!!!!!!!!!!!!!!!!')
    pause
end

fprintf('Indíce de Modulação Global: %.3f%%\n', 100*sqrt(P(Pt-P0))/P0);
fprintf('ER = %f [dB]\n', 10*log10(max(Pt)/min(Pt)))

%% Calcula Modelos
%% Modelo teórico
disp('Modelo Teórico')
TeoricoNumerico

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modelo Numérico (aproximações)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AproxNumerico

%%
CalcLogLinearEq12NewApproach

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generalização N-Tons: Aproximação linear do logaritmo Eq. Fechada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CalcEq12MatClosed

Fi = F/deltaf;
analitico2.IMPdet = 0;
ind = 1:Nc;
for k = 1:Nc % componente de produtos de intermodulação (propaga todos exceto o de interesse)
    IMP = CalcEq12RecCpp(F(ind ~= k), Fi(ind ~= k), mIM(ind ~= k), mFM(ind ~= k), dPSI(ind ~= k), phiIM(ind ~= k), beta2L, P0, N);
    IMP(abs(f) ~= F(k)) = 0;
    analitico2.IMPdet = analitico2.IMPdet + IMP;    
end
tic
Ideteqf = CalcEq12RecCpp(F, Fi, mIM, mFM, dPSI, phiIM, beta2L, P0, N);
toc

figure
fi = corta(f, indices)/1e9;
plot(fi, corta(abs(teorico.IMPdet),indices), 'sk', fi, corta(abs(aprox1.IMPdet),indices), 'or', fi, corta(abs(aprox2.IMPdet), indices), 'dm',...
    fi, corta(abs(analitico.IMPdet), indices), '*b', fi, corta(abs(analitico2.IMPdet),indices), 'pg', fi, corta(abs(loglinear.IMPdet),indices), 'sg');
legend('teórico', 'aprox1', 'aprox2', 'analitico Eq fechada', 'analitoco num', 'New approach')
axis([0 2*F(end)/1e9 0 0.5e-3])

