%% 1. TeoricoNumerico:
%		a. Calcula corrente detectada numericamente sem utilizar nenhuma
%		aproximação.
%		b. teorico.I contém todas as componentes incluindo sinal e IMPs
%       c. teorico.IMP contém apenas IMPs, porém é calculado de forma 
%       aproximada (eliminando-se um tom de cada vez e fazendo-se o cálculo).
% disp('Modelo Teórico Numérico')

teorico.At = sqrt(Pt);
teorico.dPHI = alfa/2*(log(Pt) + Adiab); 
teorico.E0 = teorico.At.*exp(1j*teorico.dPHI);
teorico.Ez = propaga(teorico.E0, L, Dhat);
teorico.I = abs(teorico.Ez).^2;
teorico.I = dc2ac(teorico.I);
teorico.Idet = t2f(teorico.I);

teorico.IMP = 0;
for k = 1:Nc % componente de produtos de intermodulação (propaga todos exceto o de interesse)
    P1 = Pt - P0*mIM(k)*cos(W(k)*t + phiIM(k));
    Adiab1 = Adiab - kappa*P0*mIM(k)/W(k)*sin(W(k)*t+phiIM(k));
    E0 = sqrt(P1).*exp(1j*alfa/2*(log(P1) + Adiab1));
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    teorico.IMP = teorico.IMP + f2t(I);    
end

teorico.IMPdet = t2f(teorico.IMP);

clear E0 Ez I P1 Adiab1