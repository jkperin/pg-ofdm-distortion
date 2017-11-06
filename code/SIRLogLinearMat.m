%% Calcula modelo analito por combinações para Eq. Orignal do artigo da Eva 
%% Peral. Este código não considera que a condição 2*f1 > 2fN é verdadeira.
%% Portanto produtos de intermodulação de qualquer ordem podem acontecer.
disp('Generalização N-Tons: Aproximação linear do logaritmo (Eva Eq. 12)');

cc3 =  [-1 -1 +1;
        -1 +1 +1;
        +1 -1 +1;
        +1 +1 -1;
        +1 +1 +1];

cc4 =  [-1 -1 -1 +1;
        -1 -1 +1 +1;
        -1 +1 -1 +1;
        -1 +1 +1 +1;
        +1 -1 -1 +1;
        +1 -1 +1 +1;
        +1 +1 -1 +1;
        +1 +1 +1 -1;
        +1 +1 +1 +1];

cc5 =  [-1 -1 -1 -1 +1;
        -1 -1 -1 +1 +1;
        -1 -1 +1 -1 +1;
        -1 -1 +1 +1 -1;
        -1 -1 +1 +1 +1;
        -1 +1 -1 -1 +1;
        -1 +1 -1 +1 +1;
        -1 +1 +1 -1 +1;
        -1 +1 +1 +1 -1;
        -1 +1 +1 +1 +1;
        +1 -1 -1 -1 +1;
        +1 -1 -1 +1 +1;
        +1 -1 +1 -1 +1;
        +1 -1 +1 +1 -1;
        +1 -1 +1 +1 +1;
        +1 +1 -1 -1 +1;
        +1 +1 -1 +1 +1;
        +1 +1 +1 -1 -1;
        +1 +1 +1 -1 +1;
        +1 +1 +1 +1 -1;
        +1 +1 +1 +1 +1];

Fi = F/deltaf;
wk = 2*pi*F;
Idet1 = zeros(1,N/2);
Idet2 = zeros(1,N/2);
Idet3 = zeros(1,N/2);
Idet4 = zeros(1,N/2);
Idet5 = zeros(1,N/2);
nk = zeros(1,Nc);
% condições fk != F[k1] && fk != F[k2] são necessárias para poder comparar
% cálculo numérico. Caso contrátio, os IMPs que caem na mesma frequência que
% os seus tons geradores serão computados este tipo de IMPs não é calculada
% no modelo numérico.

% somente calcular IMPs que caem no sinal
% 

for k1 = 1:Nc
    nk(k1) = 1;
    wimp = 2*pi*F(k1);
    thetak = -0.5*beta2L*wk*wimp;
    uk = 2*mFM.*sin(thetak);
    nn = Fi(k1) + 1;
    Idet1(nn) = Idet1(nn) + abs(P0*prod(J(nk,uk))*(1 - 1j/2*sum(mIM.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI) - J(nk+1,uk).*exp(1j*dPSI))./J(nk,uk))))^2;
    nk(k1) = 0;
    for k2 = k1+1:Nc
        fk = -F(k1)+F(k2);
        if fk >= F(1) && fk <= F(end) && fk ~= F(k1) && fk ~= F(k2)
            nn = -Fi(k1) + Fi(k2) + 1;
            wimp = 2*pi*fk;
            pok = (Fi+1 ~= nn);
            thetak = -1/2*beta2*L*wimp*wk(pok);
            uk = 2*mFM(pok).*sin(thetak);
            nk([k1,k2]) = [-1 1];
            nk2 = nk(pok);
            Idet2(nn) = Idet2(nn) + abs(P0*prod(J(nk2,uk))*(1 - 1j/2*sum(mIM(pok).*cos(thetak).*(J(nk2-1,uk).*exp(-1j*dPSI(pok)) - J(nk2+1,uk).*exp(1j*dPSI(pok)))./J(nk2,uk))))^2;
            nk([k1,k2]) = [0 0];
        end
        fk = F(k1)+F(k2);
        if fk >= F(1) && fk <= F(end) && fk ~= F(k1) && fk ~= F(k2)
            nn = Fi(k1) + Fi(k2) + 1;
            wimp = 2*pi*fk;
            pok = (Fi+1 ~= nn);
            thetak = -1/2*beta2*L*wimp*wk(pok);
            uk = 2*mFM(pok).*sin(thetak);
            nk([k1,k2]) = [1 1];
            nk2 = nk(pok);
            Idet2(nn) = Idet2(nn) + abs(P0*prod(J(nk2,uk))*(1 - 1j/2*sum(mIM(pok).*cos(thetak).*(J(nk2-1,uk).*exp(-1j*dPSI(pok)) - J(nk2+1,uk).*exp(1j*dPSI(pok)))./J(nk2,uk))))^2;
            nk([k1,k2]) = [0 0];
        end
        for k3 = k2+1:Nc
            for ii = 1:size(cc3,1)
                fk = cc3(ii,1)*F(k1) + cc3(ii,2)*F(k2) + cc3(ii,3)*F(k3);
                if fk >= F(1) && fk <= F(end) && fk ~= F(k1) && fk ~= F(k2) && fk ~= F(k3)
                    nn = cc3(ii,1)*Fi(k1) + cc3(ii,2)*Fi(k2) + cc3(ii,3)*Fi(k3) + 1;
                    wimp = 2*pi*fk;
                    pok = (Fi+1 ~= nn);
                    thetak = -1/2*beta2*L*wimp*wk(pok);
                    uk = 2*mFM(pok).*sin(thetak);
                    nk([k1,k2,k3]) = cc3(ii,:); 
                    nk2 = nk(pok);
                    Idet3(nn) = Idet3(nn) + abs(P0*prod(J(nk2,uk))*(1 - 1j/2*sum(mIM(pok).*cos(thetak).*(J(nk2-1,uk).*exp(-1j*dPSI(pok)) - J(nk2+1,uk).*exp(1j*dPSI(pok)))./J(nk2,uk))))^2;
                    nk([k1,k2,k3]) = [0 0 0];
                end
            end
        end
    end
end

% Idet(2:N/2) = conj(fliplr(Idet(N/2+2:end)));
Idet1 = [0 conj(fliplr(Idet1(2:end))) Idet1];
Idet2 = [0 conj(fliplr(Idet2(2:end))) Idet2];
Idet3 = [0 conj(fliplr(Idet3(2:end))) Idet3];
Idet4 = [0 conj(fliplr(Idet4(2:end))) Idet4];
Idet5 = [0 conj(fliplr(Idet5(2:end))) Idet5];

loglinearm.Is = f2t(Idet1);
loglinearm.Idet = Idet1;
loglinearm.IMPdet = Idet2 + Idet3 + Idet4 + Idet5;
loglinearm.IMP = f2t(loglinearm.IMPdet);

loglinearm.SIR = corta(10*log10(abs(loglinearm.Idet)./abs(loglinearm.IMPdet)), indF);

clear Idet1 Idet2 Idet3 Idet4 Idet5 It1 It2 It3 It4 It5 nk k1 k2 k3 uk thetak wimp wk Fi