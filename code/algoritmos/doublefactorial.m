function N = doublefactorial(n)
if n < -1
    disp('opção inválida')
    N = 1;
elseif n == -1 || n == 0 || n == 1
    N = 1;
else
    N = n*doublefactorial(n-2);
end