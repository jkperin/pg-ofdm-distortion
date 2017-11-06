function Idet = CalcEq12SimplOrd1Mat(F, mIM, mFM, phiIM, dPSI, theta, P0, ff, MAXIT)

J = @(n,x) besselj(n,x);

Idet = zeros(size(ff));
for nk = 1:MAXIT
        nn = find(ff == nk*F);
        u = 2*mFM.*sin(nk*theta);
        Idet(nn) = Idet(nn) + P0*1j^(nk)*exp(1j*nk*dPSI)*(J(nk,u)...
             - 1j*mIM/2.*cos(nk*theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI)))*exp(1j*nk*phiIM);
end

Idet(1:MAXIT) = fliplr(conj(Idet(MAXIT+2:end)));