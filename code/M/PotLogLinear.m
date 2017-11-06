disp('Generalização N-Tons: Cálculo com simplificações');
% Primeira ordem
Idet1c = zeros(size(f));
for k = 1:Nc % componente de sinal
    nk = zeros(size(F));
    nk(k) = 1;
    w = 2*pi*F(k);
    theta = -0.5*beta2L*2*pi*F*w;
    u = 2*mFM.*sin(theta);
    Idet1c(f == F(k)) = Idet1c(f == F(k)) + abs(P0*prod(J(nk,u))...
        *(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u))))^2;
end
Idet1c(2:N/2) = conj(fliplr(Idet1c(N/2+2:end)));
loglinear.Is = f2t(Idet1c);

fprintf('Ordem = %d\nVerificar se código C++ está calculando\n', ordem)
switch(Tipo)
    case 'normal'
        [Idet2c, Idet3c, Idet4c, Idet5c] = CalcPotLogLinear(F, F/deltaf, mIM, mFM, dPSI, beta2L, P0, N);
    case 'CAD'
        [Idet2c, Idet3c, Idet4c, Idet5c] = CalcPotLogLinearAdiabDom(F, F/deltaf, mIM, mFM, dPSI, beta2L, P0, N);
    case 'CTD'
        [Idet2c, Idet3c, Idet4c, Idet5c] = CalcPotLogLinearTransDom(F, F/deltaf, mIM, mFM, dPSI, beta2L, P0, N);
    otherwise
        disp('Opção Inválida')
end

loglinear.I2 = f2t(Idet2c);
loglinear.I3 = f2t(Idet3c);
loglinear.I4 = f2t(Idet4c);
loglinear.I5 = f2t(Idet5c);

loglinear.I = loglinear.Is;
loglinear.Idet = Idet1c;
loglinear.IMP = loglinear.I2 + loglinear.I3 + loglinear.I4 + loglinear.I5;
loglinear.IMPdet = Idet2c + Idet3c + Idet4c + Idet5c;

loglinear.Pot = corta(10*log10((loglinear.Idet + loglinear.IMPdet)/1e-3),indF);

clear Idet1c Idet2c Idet3c Idet4c Idet5c 