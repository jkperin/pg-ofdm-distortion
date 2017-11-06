function [F2, mIM2, mFM2, phiIM2, phiFM2, dPSI2, OrigInd, G] = CalcNovosTons(F, mIM, mFM, phiIM, phiFM, alfa)

Nc = length(F);
Nc2 = 2*Nc + 2*nchoosek(Nc,2);

G = 1;% - sum(mIM.^2)/16;

Fdif = zeros(1, nchoosek(Nc,2));
Fsum = zeros(size(Fdif));
phiIMdif = zeros(size(Fdif));
phiIMsum = zeros(size(Fdif));
mIMprod = zeros(size(Fdif));
mFMprod = zeros(size(Fdif));

nn = 1;
for k = 1:Nc-1
    for kk = k+1:Nc
        Fdif(nn) = F(kk)-F(k);
        Fsum(nn) = F(k)+F(kk);
        
        phiIMdif(nn) = phiIM(kk)-phiIM(k);
        phiIMsum(nn) = phiIM(k)+phiIM(kk);
         
        mIMprod(nn) = -mIM(k)*mIM(kk)/4;
        mFMprod(nn) = -alfa*mIM(k)*mIM(kk)/4; % -alfa, pois aqui alfa é definido como positivo.
        
        nn = nn + 1;
    end
end

F2 = [Fdif F Fsum 2*F];
phiIM2 = [0*phiIMdif, phiIM, 0*phiIMsum, 0*2*phiIM]; %[phiIMdif, phiIM, phiIMsum, 2*phiIM];
phiFM2 = [phiIMdif+pi/2, phiFM, phiIMsum+pi/2, 2*phiIM+pi/2];
mFM2 = [mFMprod, mFM, mFMprod, -alfa*mIM.^2/8]; % -alfa, pois aqui alfa é definido como positivo.
mIM2 = [0*mIMprod, mIM, 0*mIMprod, -0*mIM.^2/8];
dPSI2 = phiFM2 - phiIM2;
OrigInd = [zeros(size(Fdif)) ones(size(F)) zeros(size(Fsum)) zeros(size(F))];

[F2, ik] = sort(F2);
phiIM2 = phiIM2(ik);
phiFM2 = phiFM2(ik);
mFM2 = mFM2(ik);
mIM2 = mIM2(ik);
dPSI2 = dPSI2(ik);
OrigInd = OrigInd(ik);

if length(F2) ~= Nc2
    fprintf('Erro: Número de novas frequências está errado\n')
    pause
end