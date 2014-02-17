function [A] = gen_dense(eigenvalues)
    n = length(eigenvalues)
    A = zeros(n, n)
    for i = 1 : n
        c = -2.0 / n
        h = c * ones(n, 1)
        h(i) = h(i) + 1.0
        A = A + eigenvalues(i) * h * h'
    end
endfunction


function [A] = gen_tridiag(n, d, s)
    A = zeros(n, n)
    for i = 1 : n
        A(i, i) = d
    end
    for i = 1 : n - 1
        A(i, i + 1) = s
        A(i + 1, i) = s
    end
endfunction
