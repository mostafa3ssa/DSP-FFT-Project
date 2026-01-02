function X = direct_dft(x)
    N = length(x);
    
    X = zeros(N, 1);
    
    for k = 0:N-1
        for n = 0:N-1
            exponent = exp(-1j * 2 * pi * k * n / N);
            
            X(k+1) = X(k+1) + x(n+1) * exponent;
        end
    end
end