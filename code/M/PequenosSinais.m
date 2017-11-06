%%	5. PequenosSinais
%		a. Calcula corrente detectada numericamente sem utilizando aproximações de pequenos sinais
%		b. Raiz é linearizada.
%		c. Log é linearizado
%		d. Exponencial do campo é calculada com aproximação de segunda ordem
%		e. I contém todas as componentes incluindo sinal e IMPs
%		f. IMP contém apenas IMPs, porém é calculado de forma aproximada (eliminando-se um tom de cada vez e fazendo-se o cálculo).
%

smallsig.At = sqrt(P0)*(1 + AM); %sqrt(P0)*(1 + mIM/2*cos(W*t + phiIM));
smallsig.dPHI = FM; 
smallsig.E0 = smallsig.At.*(1 + 1j*smallsig.dPHI - 0.5*(smallsig.dPHI.^2));
smallsig.Ez = propaga(smallsig.E0, L, Dhat);
smallsig.I = abs(smallsig.Ez).^2;
smallsig.I = dc2ac(smallsig.I);
smallsig.Idet = t2f(smallsig.I);

smallsig.IMP = 0;
for k = 1:Nc % componente de produtos de intermodulação (propaga todos exceto o de interesse)
    At = sqrt(P0)*(1 + AM - mIM(k)/2*cos(W(k)*t + phiIM(k)));
    dPHI = FM - mFM(k)*sin(W(k)*t + phiFM(k));
    E0 = At.*(1 + 1j*dPHI - 0.5*(dPHI.^2));
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    smallsig.IMP = smallsig.IMP + f2t(I);    
end

smallsig.IMPdet = t2f(smallsig.IMP);

clear At k E0 Ez I