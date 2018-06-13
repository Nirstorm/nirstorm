function signal= Canonical(t,peakTime,uShootTime,peakDisp,uShootDisp,ratio) 
% Canonical :return the Canonical Hrf
    
    assert( isvector(t)  )

    %signal=zeros( size(t_vect) );  
    if nargin <  2, peakTime    = 6;end
    if nargin <  3, uShootTime  = 16;end
    if nargin <  4, peakDisp    = 1;end
    if nargin <  5, uShootDisp  = 1;end
    if nargin <  6, ratio       = 1/6;end
        
   signal = peakDisp^peakTime*t.^(peakTime-1).*exp(-peakDisp*t)/gamma(peakTime) - ratio*uShootDisp^uShootTime*t.^(uShootTime-1).*exp(-uShootDisp*t)/gamma(uShootTime);
   signal = signal / sum(signal);
end

