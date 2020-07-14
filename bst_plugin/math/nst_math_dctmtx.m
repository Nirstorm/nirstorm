function c = nst_math_dctmtx(n) 
    % Create a DCT matrix, see 
    % https://www.mathworks.com/help/images/ref/dctmtx.html 
    % Simimar output as: %cmat = dctmtx(nsamples)';
    
    [cc,rr] = meshgrid(0:n-1); 
    c = sqrt(2 / n) * cos(pi * (2*cc + 1) .* rr / (2 * n)); 
    c(1,:) = c(1,:) / sqrt(2); 
