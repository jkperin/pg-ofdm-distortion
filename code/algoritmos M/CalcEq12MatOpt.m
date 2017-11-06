function [Idet, I] = CalcEq12MatOpt(nkk, F, Fk, mIM, mFM, phiIM, dPSI, beta2L, P0, t)
J = @(n,x) besselj(n,x);
Wk = 2*pi*Fk;
W = 2*pi*F;

I = 0;
for k = 1:size(nkk,1)
    nk = nkk(k,:);
    nn = find(F == sum(nk.*Fk));
    theta = -(1/2)*beta2L*Wk*sum(nk.*Wk);
    u = 2*mFM.*sin(theta);
    Idet = P0*1j^(sum(nk))*exp(1j*sum(nk.*dPSI))*prod(J(nk,u))...
        *(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u)))*exp(1j*sum(nk.*phiIM));
    I = I + Idet*exp(1j*W(nn)*t) + conj(Idet)*exp(-1j*W(nn)*t);  
end

% figure(3000)
% plot(ff,THETA, 'o')
% hold on