function signal= Gamma(t,peakTime,peakDisp) 
% Gamma : apply the gamma function over t_vect
    
    %signal=zeros( size(t_vect) );  
    if nargin <  3, peakDisp = 1; end
    if nargin <  2, peakTime = 6; end

    signal = peakDisp^peakTime*t.^(peakTime-1).*exp(-peakDisp*t)/gamma(peakTime);
    signal = signal / sum(signal);
end

