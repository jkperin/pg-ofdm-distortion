%% 	3. CalcEq12MatClosed
%		a. Calcula IMPs utilizando Eq. 12 do artigo da E. Peral (Equação fechada)
%		b. Suporta somente 7 tons
%		c. Calcula as componentes de sinal e IMPs tanto em potência como em amplitude
%		d. Funções de Bessel podem ser tanto calculadas exatamente como aproximadas.
%		e. Nmax é o número máximo que os índices da equação podem atingir -Nmax <= nk <= Nmax. Nmax = 1 é para compatibilidade com o novo método de cálculo.
%		f. Esta função é utilizada para validar os códigos feitos em C.
%       g. Calcula apenas componentes na banda do sinal        
%

disp('Generalização N-Tons: Eq. Fechada');
tic
Nmax = 1;

analitico.Idet = zeros(size(F));
analitico.IMPdet = zeros(size(F));
analitico.Isdet = zeros(size(F));
analitico.IdetN = zeros(size(F));

analitico.SIR.Idet = zeros(size(F));
analitico.SIR.IMPdet = zeros(size(F));
analitico.SIR.Isdet = zeros(size(F));
analitico.SIR.IdetN = zeros(size(F));

analitico.Pot.Idet = zeros(size(F));
analitico.Pot.IMPdet = zeros(size(F));
analitico.Pot.Isdet = zeros(size(F));
analitico.Pot.IdetN = zeros(size(F));

for n1 = -Nmax:Nmax
    for n2 = -Nmax:Nmax
        for n3 = -Nmax:Nmax
            for n4 = -Nmax:Nmax
                for n5 = -Nmax:Nmax
                    for n6 = -Nmax:Nmax
                        for n7 = -Nmax:Nmax
                            nk = [n1 n2 n3 n4 n5 n6 n7];
                            fk = sum(nk.*F);
                            nn = find(fk == F);
                            if isempty(nn)
                                continue;
                            else
                                theta = -0.5*beta2L*W*sum(nk.*W);
                                u = 2*mFM.*sin(theta);
                                
                                % todas componentes amplitude
                                analitico.Idet(nn) = analitico.Idet(nn) + P0*1j^(sum(nk))*exp(1j*sum(nk.*dPSI))*prod(J(nk,u))...
                                        *(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u)))*exp(1j*sum(nk.*phiIM));
                                
                                % Todas componentes para calc da potência
                                % para tornar compatível com código c++, as
                                % IMPs devem ser calculadas até uma certa
                                % ordem
                                analitico.Pot.Idet(nn) = analitico.Pot.Idet(nn) + abs(P0*prod(J(nk,u))*(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u))))^2;
                                                                            
                                % Calc IMPs para cálculo de Pot (não compatível com modelo numérico)
                                if (length(find(nk == 0)) == 7-3 || length(find(nk == 0)) == 7-2)  % IMPs de ordem 2 e 3
                                    analitico.Pot.IMPdet(nn) = analitico.Pot.IMPdet(nn) + abs(P0*prod(J(nk,u))*(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u))))^2;
                                end
                                
                                % Calc Componente de sinal para cálculo de Pot (não compatível com modelo numérico)
                                if length(find(nk == 0)) == 7-1 % componente de sinal
                                    analitico.Pot.Isdet(nn) = analitico.Pot.Isdet(nn) + abs(P0*prod(J(nk,u))*(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u))))^2;
                                    analitico.SIR.Isdet(nn) = analitico.SIR.Isdet(nn) + abs(P0*prod(J(nk,u))*(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u))))^2;
                                end
                                
                                % Calc IMPs para cálculo de SIR (compatível com modelo numérico)
                                if (length(find(nk == 0)) == 7-3 || length(find(nk == 0)) == 7-2) && all(fk ~= abs(nk.*F)) % IMPs de ordem 2 e 3
                                    pok = 1:7;
                                    pok(fk == F) = [];
                                    nk2 = nk(pok);
                                    u2 = u(pok);
                                    mIM2 = mIM(pok);
                                    theta2 = theta(pok);
                                    dPSI2 = dPSI(pok);
                                    analitico.SIR.IMPdet(nn) = analitico.SIR.IMPdet(nn) + abs(P0*prod(J(nk2,u2))*(1 - sum(1j*mIM2/2.*cos(theta2).*(J(nk2-1, u2).*exp(-1j*dPSI2) - J(nk2+1, u2).*exp(1j*dPSI2))./J(nk2,u2)))).^2;
                                end
                                    
                                if (length(find(nk == 0)) == 4 || length(find(nk == 0)) == 5) && nk(nn) == 0
                                    analitico.IMPdet(nn) = analitico.IMPdet(nn) + P0*1j^(sum(nk(fk ~= F)))*exp(1j*sum(nk(fk ~= F).*dPSI(fk ~= F)))*prod(J(nk(fk ~= F),u(fk ~= F)))...
                                        *(1 - sum(1j*mIM(fk ~= F)/2.*cos(theta(fk ~= F)).*(J(nk(fk ~= F)-1, u(fk ~= F)).*exp(-1j*dPSI(fk ~= F)) - J(nk(fk ~= F)+1, u(fk ~= F)).*exp(1j*dPSI(fk ~= F)))./J(nk(fk ~= F),u(fk ~= F))))*exp(1j*sum(nk(fk ~= F).*phiIM(fk ~= F)));
                                elseif length(find(nk == 0)) == 6
                                    analitico.Isdet(nn) = analitico.Isdet(nn) + P0*1j^(sum(nk))*exp(1j*sum(nk.*dPSI))*prod(J(nk,u))...
                                        *(1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSI) - J(nk+1, u).*exp(1j*dPSI))./J(nk,u)))*exp(1j*sum(nk.*phiIM));
                                end
                            end
                        end
                    end
                end
            end                         
        end
    end
end

% analitico.Isdet.'
analitico.I = 0;
analitico.IMP = 0;
analitico.Is = 0;
for k = 1:length(F)
    analitico.I = analitico.I + analitico.Idet(k)*exp(1j*2*pi*F(k)*t) + conj(analitico.Idet(k))*exp(-1j*2*pi*F(k)*t);
    analitico.IMP = analitico.IMP + analitico.IMPdet(k)*exp(1j*2*pi*F(k)*t) + conj(analitico.IMPdet(k))*exp(-1j*2*pi*F(k)*t);
    analitico.Is = analitico.Is + analitico.Isdet(k)*exp(1j*2*pi*F(k)*t) + conj(analitico.Isdet(k))*exp(-1j*2*pi*F(k)*t);
end

analitico.Idet = t2f(analitico.I);
analitico.IMPdet = t2f(analitico.IMP);
analitico.Isdet = t2f(analitico.Is);

analitico.SIR.SIR = 10*log10(abs(analitico.SIR.Isdet)./abs(analitico.SIR.IMPdet));
analitico.Pot.Pot = 10*log10(analitico.Pot.Idet/1e-3);
toc