%%	2. AproxNumerico
%		a. Calcula corrente detectada numericamente sem utilizando diferentes aproximações para o logaritmo do chirp transiente
%		b. Raiz é calculada de forma exata para evitar erros inseridos pelos termos adicionais.
%		c. Aprox1. contém resultados do cálculo aproximando o log por uma aproximação linear.
%		d. Aprox2. contém resultados do cálculo aproximando o log por uma aproximação de segunda ordem
%		e. I contém todas as componentes incluindo sinal e IMPs
%		f. IMP contém apenas IMPs, porém é calculado de forma aproximada (eliminando-se um tom de cada vez e fazendo-se o cálculo).

%% Log Linear
% disp('Modelo Numérico - Log Linear (Raiz exata)')
aprox1.At = sqrt(Pt);
aprox1.dPHI = FM;
aprox1.E0 = aprox1.At.*exp(1j*aprox1.dPHI);
aprox1.Pt = abs(aprox1.E0).^2;
aprox1.Ez = propaga(aprox1.E0, L, Dhat);
aprox1.I = abs(aprox1.Ez).^2;
aprox1.I = dc2ac(aprox1.I);
aprox1.Idet = t2f(aprox1.I);

aprox1.IMP = 0;
for k = 1:Nc % componente de sinal
    At = aprox1.At - sqrt(P0)*mIM(k)/2*cos(W(k)*t + phiIM(k)); 
    dPHI = aprox1.dPHI - mFM(k)*sin(W(k)*t + phiFM(k));
    E0 = At.*exp(1j*dPHI);
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    aprox1.IMP = aprox1.IMP + f2t(I);    
end

aprox1.IMPdet = t2f(aprox1.IMP);

clear E0 At dPHI E0 Ez I

%% Log Quad
% disp('Modelo Numérico - Log Quadrático (Raiz exata)')
aprox2.At = sqrt(P0)*(1 + IM/2 - IM.^2/8);
aprox2.dPHI = alfa/2*(IM - 0*1/2*(IM).^2 + Adiab);
aprox2.E0 = aprox2.At.*exp(1j*aprox2.dPHI);
aprox2.Ez = propaga(aprox2.E0, L, Dhat);
aprox2.I = abs(aprox2.Ez).^2;
aprox2.I = dc2ac(aprox2.I);
aprox2.Idet = t2f(aprox2.I);

aprox2.IMP = 0;
for k = 1:Nc % componente de sinal
    IM1 = IM - mIM(k)*cos(W(k)*t + phiIM(k));
    Adiab1 = Adiab - kappa*P0*mIM(k)/W(k)*sin(W(k)*t+phiIM(k));
    dPHI = alfa/2*(IM1 - 1/2*(IM1).^2 + Adiab1 + log(P0));
    At = sqrt(P0)*(1 + IM1/2 - IM1.^2/8);
    E0 = At.*exp(1j*dPHI);
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    aprox2.IMP = aprox2.IMP + f2t(I);    
end

aprox2.IMPdet = t2f(aprox2.IMP);

clear E0 At dPHI E0 Ez I
