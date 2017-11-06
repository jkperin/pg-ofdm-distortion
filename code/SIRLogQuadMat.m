%% Calcula modelo analito por combinações para teoria de grandes sinais com
%% aproximação de segunda ordem para o logaritmo do chirp transiente.
%% Este código não considera que a condição 2*f1 > 2fN é verdadeira.
%% Portanto produtos de intermodulação de qualquer ordem podem acontecer.
disp('Generalização N-Tons: Aproximação quadrática do logaritmo');

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
    
Nc2 = 2*(Nc) + 2*nchoosek(Nc,2);
[F2, mIM2, mFM2, phiIM2, phiFM2, dPSI2, OrigInd, G] = CalcNovosTons(F, mIM, mFM, phiIM, phiFM, alfa);

%%
phi2 = 0;
A2 = G;
for k = 1:Nc2
    phi2 = phi2 + mFM2(k)*sin(2*pi*F2(k)*t + phiFM2(k));
    A2 = A2 + mIM2(k)/2*cos(2*pi*F2(k)*t + phiIM2(k));
end

phi = 0;
A1 = 1;
for k = 1:Nc
    phi = phi + mFM(k)*sin(2*pi*F(k)*t + phiFM(k));
    A1 = A1 + mIM(k)/2*cos(2*pi*F(k)*t + phiIM(k));
end


phiteo = alfa/2*(log(1 + IM) + Adiab);
phi1teo = alfa/2*(IM + Adiab);
phi2teo = alfa/2*(IM - IM.^2/2 + Adiab);

Ateo = sqrt(1 + IM);
A1teo = 1 + IM/2;
A2teo = 1 + IM/2 - IM.^2/8;

figure
plot(t, dc2ac(phi2teo), t, dc2ac(phi), t, dc2ac(phi2),'--')
legend('teórico', 'linear', 'quad')
title('Fase')
figure
plot(t, A2teo, t, A1, t, A2,'--') 
legend('teórico', 'linear', 'quad')
title('Amplitude')
%%
Idet1a = zeros(1,N/2);
Idet2a = zeros(1,N/2);
Idet3a = zeros(1,N/2);
Idet4a = zeros(1,N/2);
Idet5a = zeros(1,N/2);

Fi = F2/deltaf;
wk = 2*pi*F2;
nk = zeros(size(F2));

for k1 = 1:Nc2
    if F2(k1) >= F(1) && F2(k1) <= F(end)
        nn = Fi(k1) + 1;
        wimp = 2*pi*F2(k1);
        thetak = -1/2*beta2*L*wimp*wk;
        uk = 2*mFM2.*sin(thetak);
        nk(k1) = 1;
        Idet1a(nn) = Idet1a(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
        nk(k1) = 0;
    end
    for k2 = k1+1:Nc2
        fk = -F2(k1)+F2(k2);
        if fk >= F(1) && fk <= F(end)
            nn = -Fi(k1) + Fi(k2) + 1;
            wimp = 2*pi*fk;
            thetak = -1/2*beta2*L*wimp*wk;
            uk = 2*mFM2.*sin(thetak);
            nk([k1,k2]) = [-1 1];
            Idet2a(nn) = Idet2a(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
            nk([k1,k2]) = [0 0];
        end
        fk = F2(k1)+F2(k2);
        if fk >= F(1) && fk <= F(end)
            nn = Fi(k1) + Fi(k2) + 1;
            wimp = 2*pi*fk;
            pok = (Fi+1 ~= nn);
            thetak = -1/2*beta2*L*wimp*wk;
            uk = 2*mFM2.*sin(thetak);
            nk([k1,k2]) = [1 1];
            Idet2a(nn) = Idet2a(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
            nk([k1,k2]) = [0 0];
        end
        for k3 = k2+1:Nc2
            for ii = 1:size(cc3,1)
                fk = cc3(ii,1)*F2(k1) + cc3(ii,2)*F2(k2) + cc3(ii,3)*F2(k3);
                if fk >= F(1) && fk <= F(end)
                    nn = cc3(ii,1)*Fi(k1) + cc3(ii,2)*Fi(k2) + cc3(ii,3)*Fi(k3) + 1;
                    wimp = 2*pi*fk;
                    thetak = -1/2*beta2*L*wimp*wk;
                    uk = 2*mFM2.*sin(thetak);
                    nk([k1,k2,k3]) = cc3(ii,:); 
                    Idet3a(nn) = Idet3a(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
                    nk([k1,k2,k3]) = [0 0 0];
                end
            end
%             for k4 = k3+1:Nc2
%                 for ii = 1:size(cc4,1)
%                     fk = cc4(ii,1)*F2(k1) + cc4(ii,2)*F2(k2) + cc4(ii,3)*F2(k3) + cc4(ii,4)*F2(k4);
%                     if fk >= F(1) && fk <= F(end) && fk ~= F2(k1) && fk ~= F2(k2) && fk ~= F2(k3) && fk ~= F2(k4)
%                         nn = cc4(ii,1)*Fi(k1) + cc4(ii,2)*Fi(k2) + cc4(ii,3)*Fi(k3) + cc4(ii,4)*Fi(k4) +  1;
%                         wimp = 2*pi*fk;
%                         pok = (Fi+1 ~= nn);
%                         thetak = -1/2*beta2*L*wimp*wk;
%                         uk = 2*mFM2.*sin(thetak);
%                         nk([k1,k2,k3,k4]) = cc4(ii,:); 
%                         nk = nk;
%                         Idet4(nn) = Idet4(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
%                         nk([k1,k2,k3,k4]) = [0 0 0 0];
%                     end
%                 end
%             end
        end
    end
end    

Idet1a = [0 conj(fliplr(Idet1a(2:end))) Idet1a];
Idet2a = [0 conj(fliplr(Idet2a(2:end))) Idet2a];
Idet3a = [0 conj(fliplr(Idet3a(2:end))) Idet3a];
Idet4a = [0 conj(fliplr(Idet4a(2:end))) Idet4a];
Idet5a = [0 conj(fliplr(Idet5a(2:end))) Idet5a];

clear F2 mIM2 mFM2  phiIM2 phiFM2 dPSI2 Fi wk uk thetak nn wimp nk 

Idet1 = zeros(1,N/2);
Idet2 = zeros(1,N/2);
Idet3 = zeros(1,N/2);
Idet4 = zeros(1,N/2);
Idet5 = zeros(1,N/2);

Nc2 = 2*(Nc-1) + 2*nchoosek(Nc-1,2);

for kn = 1:Nc
    pok = 1:Nc; pok(kn) = [];
    [F2, mIM2, mFM2, phiIM2, phiFM2, dPSI2, OrigInd, G] = CalcNovosTons(F(pok), mIM(pok), mFM(pok), phiIM(pok), phiFM(pok), alfa);
    Fi = F2/deltaf;
    wk = 2*pi*F2;
    nk = zeros(1,Nc2);
    % condições fk != F[k1] && fk != F[k2] são necessárias para poder comparar
    % cálculo numérico. Caso contrátio, os IMPs que caem na mesma frequência que
    % os seus tons geradores serão computados este tipo de IMPs não é calculada
    % no modelo numérico.

    for k1 = 1:Nc2
        fk = F2(k1);
        if fk == F(kn)
            nn = Fi(k1) + 1;
            wimp = 2*pi*fk;
            thetak = -1/2*beta2*L*wimp*wk;
            uk = 2*mFM2.*sin(thetak);
            nk(k1) = 1;
            Idet1(nn) = Idet1(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
            nk(k1) = 0;
        end
        for k2 = k1+1:Nc2
            fk = -F2(k1)+F2(k2);
            if fk == F(kn)
                nn = -Fi(k1) + Fi(k2) + 1;
                wimp = 2*pi*fk;
                thetak = -1/2*beta2*L*wimp*wk;
                uk = 2*mFM2.*sin(thetak);
                nk([k1,k2]) = [-1 1];
                Idet2(nn) = Idet2(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
                nk([k1,k2]) = [0 0];
            end
            fk = F2(k1)+F2(k2);
            if fk == F(kn)
                nn = Fi(k1) + Fi(k2) + 1;
                wimp = 2*pi*fk;
                pok = (Fi+1 ~= nn);
                thetak = -1/2*beta2*L*wimp*wk;
                uk = 2*mFM2.*sin(thetak);
                nk([k1,k2]) = [1 1];
                Idet2(nn) = Idet2(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
                nk([k1,k2]) = [0 0];
            end
            for k3 = k2+1:Nc2
                for ii = 1:size(cc3,1)
                    fk = cc3(ii,1)*F2(k1) + cc3(ii,2)*F2(k2) + cc3(ii,3)*F2(k3);
                    if fk == F(kn)
                        nn = cc3(ii,1)*Fi(k1) + cc3(ii,2)*Fi(k2) + cc3(ii,3)*Fi(k3) + 1;
                        wimp = 2*pi*fk;
                        thetak = -1/2*beta2*L*wimp*wk;
                        uk = 2*mFM2.*sin(thetak);
                        nk([k1,k2,k3]) = cc3(ii,:); 
                        Idet3(nn) = Idet3(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
                        nk([k1,k2,k3]) = [0 0 0];
                    end
                end
%                 for k4 = k3+1:Nc2
%                     for ii = 1:size(cc4,1)
%                         fk = cc4(ii,1)*F2(k1) + cc4(ii,2)*F2(k2) + cc4(ii,3)*F2(k3) + cc4(ii,4)*F2(k4);
%                         if fk >= F(1) && fk <= F(end) && fk ~= F2(k1) && fk ~= F2(k2) && fk ~= F2(k3) && fk ~= F2(k4)
%                             nn = cc4(ii,1)*Fi(k1) + cc4(ii,2)*Fi(k2) + cc4(ii,3)*Fi(k3) + cc4(ii,4)*Fi(k4) +  1;
%                             wimp = 2*pi*fk;
%                             pok = (Fi+1 ~= nn);
%                             thetak = -1/2*beta2*L*wimp*wk;
%                             uk = 2*mFM2.*sin(thetak);
%                             nk([k1,k2,k3,k4]) = cc4(ii,:); 
%                             Idet4(nn) = Idet4(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
%                             nk([k1,k2,k3,k4]) = [0 0 0 0];
%                         end
%                     end
%                 end
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

clear It1 It2 It3 It4 It5 nk k1 k2 k3 uk thetak wimp wk Fi

%%
logquad.IMPdet = Idet1 + Idet2 + Idet3 + Idet4 + Idet5;
logquad.IMP = f2t(logquad.IMPdet);

logquad.Idet = Idet1a + Idet2a + Idet3a + Idet4a + Idet5a - logquad.IMPdet;
logquad.Is = f2t(logquad.Idet);

logquad.SIR = corta(10*log10(abs(logquad.Idet)./abs(logquad.IMPdet)), indF);


