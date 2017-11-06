%% Calcula modelo analito por combinações para teoria de grandes sinais com
%% aproximação de segunda ordem para o logaritmo do chirp transiente.
%% Este código não considera que a condição 2*f1 > 2fN é verdadeira.
%% Portanto produtos de intermodulação de qualquer ordem podem acontecer.
disp('Generalização N-Tons: Aproximação quadrática do logaritmo');
tic
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
        

Idet1sm = zeros(1, N/2);
Idet1m = zeros(1, N/2);
Idet2m = zeros(1, N/2);
Idet3m = zeros(1, N/2);
Idet4m = zeros(1, N/2);
Idet5m = zeros(1, N/2);
    
Nc2 = 2*(Nc) + 2*nchoosek(Nc,2);
[F2, mIM2, mFM2, ~, ~, dPSI2, OrigInd, G] = CalcNovosTons(F, mIM, mFM, phiIM, phiFM, alfa);

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
        if OrigInd(k1) % componente de sinal
            Idet1sm(nn) = Idet1sm(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
        else % IMP de ordem 1
            Idet1m(nn) = Idet1m(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
        end
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
            Idet2m(nn) = Idet2m(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
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
            Idet2m(nn) = Idet2m(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
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
                    Idet3m(nn) = Idet3m(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
                    nk([k1,k2,k3]) = [0 0 0];
                end
            end
%             for k4 = k3+1:Nc2
%                 for ii = 1:size(cc4,1)
%                     fk = cc4(ii,1)*F2(k1) + cc4(ii,2)*F2(k2) + cc4(ii,3)*F2(k3) + cc4(ii,4)*F2(k4);
%                     if fk >= F(1) && fk <= F(end)
%                         nn = cc4(ii,1)*Fi(k1) + cc4(ii,2)*Fi(k2) + cc4(ii,3)*Fi(k3) + cc4(ii,4)*Fi(k4) +  1;
%                         wimp = 2*pi*fk;
%                         pok = (Fi+1 ~= nn);
%                         thetak = -1/2*beta2*L*wimp*wk;
%                         uk = 2*mFM2.*sin(thetak);
%                         nk([k1,k2,k3,k4]) = cc4(ii,:); 
%                         Idet4m(nn) = Idet4m(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;
%                         nk([k1,k2,k3,k4]) = [0 0 0 0];
%                     end
%                 end
%             end
        end
    end
end    

Idet1sm = [0 conj(fliplr(Idet1sm(2:end))) Idet1sm];
Idet1m = [0 conj(fliplr(Idet1m(2:end))) Idet1m];
Idet2m = [0 conj(fliplr(Idet2m(2:end))) Idet2m];
Idet3m = [0 conj(fliplr(Idet3m(2:end))) Idet3m];
Idet4m = [0 conj(fliplr(Idet4m(2:end))) Idet4m];
Idet5m = [0 conj(fliplr(Idet5m(2:end))) Idet5m];

logquad.Idet = Idet1sm;
logquad.Is = f2t(logquad.Idet);

logquad.IMPdet = Idet1m + Idet2m + Idet3m + Idet4m + Idet5m;
logquad.IMP = f2t(logquad.IMPdet);

logquad.Pot = corta(10*log10((logquad.Idet + logquad.IMPdet)/1e-3), indF);

clear Idet1m Idet2m Idet3m Idet4m Idet5m nk k1 k2 k3 uk thetak wimp wk Fi
toc