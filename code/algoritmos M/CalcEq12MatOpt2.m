function I = CalcEq12MatOpt2(nkk, F, mIM, mFM, phiIM, dPSI, beta2L, P0, G, t)
J = @(n,x) besselj(n,x);

I = 0;
for k = 1:size(nkk,1)
    nk = nkk(k,:);
    w = 2*pi*sum(nk.*F);
    theta = -(1/2)*beta2L*2*pi*F*w;
    u = 2*mFM.*sin(theta);
    Idet = P0*1j^(sum(nk))*exp(1j*sum(nk.*dPSI))*prod(J(nk,u))...
        *G*(G - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u)))*exp(1j*sum(nk.*phiIM));
    I = I + Idet*exp(1j*w*t) + conj(Idet)*exp(-1j*w*t);  
end
