function [Idet, I] = CalcEq12SimplOrd3Mat(F, mFM, phiIM, beta2L, P0, alfa, t, Deltaf, MAXIT)

W = 2*pi*F;

ffmax = MAXIT*(sum(F));
ff = -abs(ffmax):Deltaf:abs(ffmax);
ww = 2*pi*ff;

Idet = zeros(size(ff));
for n1 = [-MAXIT:-1 1:MAXIT]
    for n2 = [-MAXIT:-1 1:MAXIT]
        for n3 = [-MAXIT:-1 1:MAXIT]
            nk = [n1 n2 n3];
            if sum(nk.*F) > 0
                nn = find(ff == sum(nk.*F));
            else
                continue;
            end
            theta = -(1/2)*beta2L*W*sum(nk.*W);
            u = 2*mFM.*sin(theta);
%             Idet(nn) = Idet(nn) + P0*prod(-u.*nk/2)*(1 - 0*1/alfa*sum(nk./tan(theta)))*exp(1j*sum(nk.*phiIM));
            Idet(nn) = Idet(nn) + P0*prod(-u.*nk/2)*exp(1j*sum(nk.*phiIM));
        end
    end
end

Idet(ff < 0) = fliplr(conj(Idet(ff > 0)));

I = 0;
for k = 1:(length(ff)-1)/2
    I = I + Idet(k)*exp(1j*ww(k)*t) + Idet(end-k+1)*exp(-1j*ww(k)*t);
end
