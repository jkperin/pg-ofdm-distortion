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
teorico.dPHI = mFM*sin(W*t + phiFM);
teorico.E0 = teorico.At.*exp(1j*teorico.dPHI);
Dhat = -1j*0.5*beta2*(2*pi*ifftshift(f)).^2;
teorico.Ez = propaga(teorico.E0, L, Dhat);
teorico.I = abs(teorico.Ez).^2;
teorico.Idet = t2f(teorico.I);

teorico.Transiente = alfa/2*log(Pt);
teorico.Adiabatico = alfa/2*x;