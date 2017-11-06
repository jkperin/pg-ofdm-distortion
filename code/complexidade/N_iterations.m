% clear, clc, close all
n1 = 0;
n2 = 0;
n3 = 0;
n4 = 0;
n5 = 0;
% Nc = 100;

%% Fórmulas analíticas
% n2 = (N-2) + (N-3) + ... + 1 = 1/2*N*(N-1)
% n3 = (N-1)*(N-2) + (N-2)*(N-3) + (N-3)*(N-4) + ... + 1 (N-2) vezes)
% = (N-4)*N^2/2 - N*(2*N+4)*(N-4)/4 + 1/3*((N-4)^2+9*(N-4)+26)*(N-4)/2
% = 1/6*N*(N-1)*(N-2)
% n4 = 1/24*N*(N-1)*(N-2)*(N-3)
% n5 = 1/120*N*(N-1)*(N-2)*(N-3)*(N-4)
% nk = \frac{1}{k!}\prod_{j=0}^{k-1} (N - 2(k-1) + j)
for k1 = 1:Nc
    n1 = n1 + 1;
    for k2 = k1+1:Nc
        n2 = n2 + 1;
        for k3 = k2+1:Nc
            n3 = n3 + 1;
%             for k4 = k3+1:Nc
%                 n4 = n4 + 1;
%                 for k5 = k4+1:Nc
%                     n5 = n5 + 1;
%                 end
%             end
        end
    end
end
n = n1 + n2*2 + n3*5 + n4*9 + n5*21
1521.793205/97735136
(n*250.002120/97735136)/3600

% function N = FWMterms(Nc)
% if mod(Nc,2) == 0
%     N = (Nc^2 -6*Nc + 2*Nc*(1:Nc) - 2*(1:Nc).^2 + 2*(1:Nc) + 4)/4;
% else
%     N = (Nc^2 -6*Nc + 2*Nc*(1:Nc) - 2*(1:Nc).^2 + 2*(1:Nc) + 4 + (-1).^(1:Nc))/4;
% end