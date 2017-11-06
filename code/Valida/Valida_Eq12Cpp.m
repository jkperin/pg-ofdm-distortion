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
Nc = 7;           % Número de subportadoras
fc = 5e9;         % Frequência central (Hz)
BW = Nc*1e9;           % Bandwidth (Hz)
Deltaf = BW/Nc;     

F = (-3:3)*Deltaf + fc;
W = 2*pi*F;

% if 2*min(F) > max(F)
%     disp('Sistema Suboitava, não há distorção harmônica');
%     pause
% end

%% laser
mIM = 0.05*ones(size(F));
QPSK = [pi/4 3*pi/4 -3*pi/4 -pi/4];
PSK8 = [0 pi/4 pi/2 3*pi/4 pi -3*pi/4 -pi/2 -pi/4];
phiIM = QPSK(randi(3, [1, Nc]) + 1);

lambda = 1550e-9;
c = 299792458;

% Artigo Tang (30Gb/s Signal transmission over 40km DML...)
laserfreq = c/lambda;
alfa = 3;   
OpticalConfinementEff = 0.07;
GainCompressionFac = 7.5e-23;
Vact = 29.2e-18;
DifQuantumEffpFacet = 0.139;
hPlanck = 6.626068e-34;
kappa = 0*2*OpticalConfinementEff*GainCompressionFac/(Vact*DifQuantumEffpFacet*hPlanck*laserfreq);

P0 = 4.3e-3;

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
L = 85.5e3;

beta2 = -D*lambda^2/(2*pi*c);
theta = -0.5*beta2*W.^2*L;

%%
Pt = P0 + zeros(size(t));
x = zeros(size(t));

for k = 1:Nc
    Pt = Pt + P0*mIM(k)*cos(W(k)*t+phiIM(k));
    x = x + kappa*P0*mIM(k)/W(k)*sin(W(k)*t+phiIM(k));
end

if min(Pt) <= 0
    disp('erro! Potência negativa!!!!!!!!!!!!!!!!!!')
    pause
end

fprintf('Indíce de Modulação Geral(+): %.3f\n', max(Pt)/mean(Pt) - 1);
fprintf('Indíce de Modulação Geral(-): %.3f\n', 1 - min(Pt)/mean(Pt));
fprintf('ER = %f [dB]\n', 10*log10(max(Pt)/min(Pt)))

MAXIT = 3;

%% Calcula Modelos
%% Modelo teórico
TeoricoNumerico

%%
Nmax = MAXIT;
ffmin = -Nmax*(sum(F));
ffmax = (Nmax-1)*(sum(F));
dff = Deltaf;
f3 = ffmin:dff:ffmax;

mFM = zeros(size(W));
phiFM = zeros(size(W));
dPSI = zeros(size(W));
for i = 1:length(W)
    FM = alfa/2*(1 + kappa*P0./(1j*W(i)))*1j*mIM(i)*exp(1j*phiIM(i));
    mFM(i) = abs(FM);
    phiFM(i) = angle(FM);
    dPSI(i) = phiFM(i) - phiIM(i);  
end

%% Matlab eq. 12
Nmax = MAXIT;
ffmin = -Nmax*(sum(F));
ffmax = (Nmax-1)*(sum(F));
dff = Deltaf;
ff = ffmin:dff:ffmax;

eva12.Idet = zeros(size(ff));
eva12.I = zeros(size(t));
tic
for n1 = -Nmax:Nmax-1
    for n2 = -Nmax:Nmax-1
        for n3 = -Nmax:Nmax-1
            for n4 = -Nmax:Nmax-1
                for n5 = -Nmax:Nmax-1
                    for n6 = -Nmax:Nmax-1
                        for n7 = -Nmax:Nmax-1      
                            nk = [n1 n2 n3 n4 n5 n6 n7];
                            nn = sum(nk.*F/dff) + sum(Nmax*F/dff) + 1;
                            theta = -(1/2)*beta2*L*W*sum(nk.*W);
                            u = 2*mFM.*sin(theta);
                            eva12.Idet(nn) = eva12.Idet(nn) + P0*1j^(sum(nk))*exp(1j*sum(nk.*dPSI))*prod(J(nk,u))...
                                *(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u)));
                        end
                    end
                end
            end                         
        end
    end
end

eva12.Idet(floor(length(eva12.Idet)/2) +1) = 0;

toc
%% C++ eq 12
tic
cpp.Idet = CalcEq12RecCpp(length(f3), mIM, mFM, dPSI, W, F/dff, beta2, L, P0, Nmax);
toc
%% C++ eq 12
tic
f4 = -ffmax:dff:ffmax-dff;
cppopt.Idet = CalcEq12RecCppOpt(length(f4), mIM, mFM, dPSI, W, F/dff, beta2, L, P0, Nmax, f4(end)/dff);
toc
%% Gráficos
plot(corta(f, indices), corta(abs(teorico.Idet), indices), 'o', ff, abs(eva12.Idet), '*m', f3, abs(cpp.Idet), 'sk', f4, abs(cppopt.Idet), 'vr')
legend('Teorico', 'Matlab', 'Cpp', 'Cpp Opt')
axis([-30e9 30e9 0 0.1*P0])