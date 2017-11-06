%% calcula nk
function nk = calc_N_it(k, N)
    produt = ones(size(N));
    for j = 0:k-1
        produt = produt.*(N - j);
    end
    nk = produt/factorial(k);
    