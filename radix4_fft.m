function X = radix4_fft(x)
    N = length(x);
    
    if N == 1
        X = x;
        return;
    end
    
    if mod(log2(N)/2, 1) ~= 0
        error('Input length N must be a power of 4 (e.g., 16, 64, 256)');
    end
    
    x_0 = x(1:4:end);
    x_1 = x(2:4:end);
    x_2 = x(3:4:end);
    x_3 = x(4:4:end);
    
    F0 = radix4_fft(x_0);
    F1 = radix4_fft(x_1);
    F2 = radix4_fft(x_2);
    F3 = radix4_fft(x_3);
    
    X = zeros(N, 1);
    M = N / 4; 
    
    k = (0 : M-1)';
    
    W1 = exp(-1j * 2 * pi * k / N);     
    W2 = exp(-1j * 2 * pi * 2 * k / N); 
    W3 = exp(-1j * 2 * pi * 3 * k / N); 
    
    T1 = W1 .* F1;
    T2 = W2 .* F2;
    T3 = W3 .* F3;
    
    X(1:M)       = F0 + T1 + T2 + T3;           % First Quarter
    X(M+1:2*M)   = F0 - 1j*T1 - T2 + 1j*T3;     % Second Quarter
    X(2*M+1:3*M) = F0 - T1 + T2 - T3;           % Third Quarter
    X(3*M+1:4*M) = F0 + 1j*T1 - T2 - 1j*T3;     % Fourth Quarter
end