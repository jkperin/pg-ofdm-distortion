clear, clc, close all

NC = 100:100:800;

for kk = 1:length(NC)
    Nc = NC(kk);
    N_iterations
    nn(kk) = n;
end

analitico = calc_N_it(1,NC) + 2*calc_N_it(2,NC) + 5*calc_N_it(3,NC);
plot(NC, nn, NC, analitico, 'o')

