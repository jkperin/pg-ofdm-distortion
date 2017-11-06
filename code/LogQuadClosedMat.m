disp('Generalização N-Tons: Log Quad Eq. Fechada');

Nc2 = 2*(Nc) + 2*nchoosek(Nc,2);
[F2, mIM2, mFM2, phiIM2, phiFM2, dPSI2, OrigInd, G] = CalcNovosTons(F, mIM, mFM, phiIM, phiFM, alfa);

if Nc2 ~= 12
    disp('Erro: Número de subportadoras não permitido')
    pause
end

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

analitico2.Idet = zeros(1,N);

Nmax = 2;

for n1 = -Nmax:Nmax
    for n2 = -Nmax:Nmax
        for n3 = -Nmax:Nmax
            for n4 = -Nmax:Nmax
                for n5 = -Nmax:Nmax
                    for n6 = -Nmax:Nmax
                        for n7 = -Nmax:Nmax
                            for n8 = -Nmax:Nmax
                                for n9 = -Nmax:Nmax
                                    for n10 = -Nmax:Nmax
                                        for n11 = -Nmax:Nmax
                                            for n12 = -Nmax:Nmax
                                                nk = [n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11 n12];
                                                fk = sum(nk.*F2);
                                                nn = sum(nk.*Fi) + N/2 + 1;
                                                thetak = -0.5*beta2L*wk*sum(nk.*wk);
                                                uk = 2*mFM2.*sin(thetak);

                                                analitico2.Idet(nn) = analitico2.Idet(nn) + abs(P0*prod(J(nk,uk))*G*(G - 1j/2*sum(mIM2.*cos(thetak).*(J(nk-1,uk).*exp(-1j*dPSI2) - J(nk+1,uk).*exp(1j*dPSI2))./J(nk,uk))))^2;                                        
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end                         
        end
    end
end

analitico2.I = f2t(analitico2.Idet);
analitico2.Pot = 10*log10(corta(analitico2.Idet,indF)/1e-3);

clear n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 Idet Ateo A1teo A2teo phiteo phi1teo phi2teo