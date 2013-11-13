function [A] = f(d)
    n = length(d)
    A = zeros(n, n)
    for i = 1 : n
        c = -2.0 / n
        h = c * ones(n, 1)
        h(i) = h(i) + 1.0
        A = A + d(i) * h * h'
    end
endfunction
