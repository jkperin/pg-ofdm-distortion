%% Calcula numericamente a corrente detectada para as aproximações linear e
%% quadrática do logaritmo. A raiz quadrada é aproximada linearmente

disp('Modelo Log Linear Numérico (Raiz Linearizada)')
aprox1.At = sqrt(P0)*(1 + AM);
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
    x1 = aprox1.dPHI - mFM(k)*sin(W(k)*t + phiFM(k));
    E0 = At.*exp(1j*x1);
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    aprox1.IMP = aprox1.IMP + f2t(I);    
end

aprox1.IMPdet = t2f(aprox1.IMP);

%% Modelo Log Quadrático Numérico (Raiz Linearizada)
disp('Modelo Log Quadrático Numérico (Raiz Linearizada)')
aprox2.At = sqrt(P0)*(1 + AM);
aprox2.dPHI = alfa/2*(IM - 1/2*(-IM).^2 + x + log(P0));
aprox2.E0 = aprox2.At.*exp(1j*aprox2.dPHI);
aprox2.Ez = propaga(aprox2.E0, L, Dhat);
aprox2.I = abs(aprox2.Ez).^2;
aprox2.I = dc2ac(aprox2.I);
aprox2.Idet = t2f(aprox2.I);

aprox2.IMP = 0;
for k = 1:Nc % componente de sinal
    At = aprox2.At - sqrt(P0)*mIM(k)/2*cos(W(k)*t + phiIM(k));
    x1 = IM - mIM(k)*cos(W(k)*t + phiIM(k));
    x1 = alfa/2*(x1 - 1/2*(-x1).^2 + log(P0));
    E0 = At.*exp(1j*x1);
    Ez = propaga(E0, L, Dhat);
    I = abs(Ez).^2;
    I = t2f(I);
    I(abs(f) ~= F(k)) = 0;
    aprox2.IMP = aprox2.IMP + f2t(I);    
end

aprox2.IMPdet = t2f(aprox2.IMP);