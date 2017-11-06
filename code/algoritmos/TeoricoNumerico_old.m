%% Calcula modelo teórico numericamente.
% E(t,0) = sqrt(Pt)*exp(1j*alfa/2*(log(Pt) + int(Pt)))

if length(W) == 1
    FM = alfa/2*(1 + kappa*P0/(1j*W))*1j*mIM*exp(1j*phiIM);
    mFM = abs(FM);
    phiFM = angle(FM);
    dPSI = phiFM - phiIM;
end

disp('Modelo Teórico Numérico')
teorico.str = 'Modelo Teórico Numérico';
teorico.At = sqrt(Pt); %sqrt(P0)*(1 + mIM/2*cos(W*t + phiIM));
teorico.dPHI = alfa/2*(log(Pt) + x); 
teorico.E0 = teorico.At.*exp(1j*teorico.dPHI);
Dhat = -1j*0.5*beta2*(2*pi*ifftshift(f)).^2;
teorico.Ez = propaga(teorico.E0, L, Dhat);
teorico.I = abs(teorico.Ez).^2;
teorico.Idet = t2f(teorico.I);

teorico.Is = 0;
for k = 1:Nc
    P1 = P0*(1 + mIM(k)*cos(W(k)*t + phiIM(k)));
    x1 = kappa*P0*mIM(k)/W(k)*sin(W(k)*t+phiIM(k));
    E0 = sqrt(P1).*exp(1j*alfa/2*(log(P1) + x1));
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    teorico.Is = teorico.Is + f2t(I);
end

teorico.IMP = 0;
for k = 1:Nc % componente de sinal
    P1 = Pt - P0*mIM(k)*cos(W(k)*t + phiIM(k));
    x1 = kappa*P0*mIM(k)/W(k)*sin(W(k)*t+phiIM(k));
    E0 = sqrt(P1).*exp(1j*alfa/2*(log(P1) + x1));
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    teorico.IMP = teorico.IMP + f2t(I);    
end

teorico.IMPdet = t2f(teorico.IMP);

teorico.Transiente = alfa/2*log(Pt);
teorico.Adiabatico = alfa/2*x;