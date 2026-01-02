function X = radix2_fft(x)
    N = length(x);
    
    if N == 1
        X = x;
        return;
    end
    
    if mod(log2(N), 1) ~= 0
        error('Input length N must be a power of 2 for Radix-2 FFT');
    end
    
    x_even = x(1:2:end);
    x_odd  = x(2:2:end);
    
    X_even = radix2_fft(x_even);
    X_odd  = radix2_fft(x_odd);
    
    X = zeros(N, 1);
    
    k = (0 : N/2 - 1)'; 
    W = exp(-1j * 2 * pi * k / N);
    
    X(1 : N/2)     = X_even + W .* X_odd;
    X(N/2 + 1 : N) = X_even - W .* X_odd;
end