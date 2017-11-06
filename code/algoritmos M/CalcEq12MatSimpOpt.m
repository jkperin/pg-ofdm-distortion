function [Idet, I] = CalcEq12MatSimpOpt(nkk, F, Fk, mFM, phiIM, beta2L, P0, t)
Wk = 2*pi*Fk;
W = 2*pi*F;

I = 0;
for k = 1:size(nkk,1)
    nk = nkk(k,:);
    nn = find(F == sum(nk.*Fk));
    theta = -(1/2)*beta2L*Wk*sum(nk.*Wk);
    u = 2*mFM.*sin(theta);
    Idet = P0*prod(-u.*nk/2)*exp(1j*sum(nk.*phiIM));%*(1 - 1/5*sum(nk./tan(theta)));
    I = I + Idet*exp(1j*W(nn)*t) + conj(Idet)*exp(-1j*W(nn)*t);  
end

% figure(3000)
% plot(ff,THETA, 'o')
% hold on