function [COV, BSL] = be_randomize_ridge_baseline( B, O, varargin )          

% Initialize variable
warning('OFF')
[nbS, nbT]  = size(B);
nbT         = nbT - mod(nbT,2);
method      = 'average_covariance';
if nargin>2
    method  = varargin{1};
end

% Loop on sensors
for ll = 1 : size(B,1)
    [bsl,O]	=   be_cwsparse(B(ll,:), O);
    FRQlims =   be_closest(O.ridges.frequency_range,O.wavelet.freqs_analyzed);
    
    if ll==1
        BSL = zeros(nbS,diff(FRQlims)+1,size(B,2));
    end
    
    BSL(ll,:,:) =   bsl( FRQlims(1):FRQlims(2),: );
end

% Compute covariance
switch method
    
    case 'none'
        
        COV = [];
        return
    
    case 'average_baseline'
    
        b = zeros( nbS, nbT/2 );
        for ii = 1 : 100
            st = fix(rand*nbT/2)+1;
            nd = st + nbT/2-1;
            sc = fix(rand*(diff(FRQlims)+1))+1;
            b  = b + squeeze( BSL(:,sc,st:nd) )/100;    
        end

        COV = cov( real(b') ) + 1i*cov( imag(b') );
    
    case 'average_covariance'

        COV = zeros( nbS );
        for ii = 1 : 100
            st = fix(rand*nbT/2)+1;
            nd = st + nbT/2-1;
            sc = fix(rand*(diff(FRQlims)+1))+1;
            b  = squeeze( BSL(:,sc,st:nd) );
            COV= COV + ( cov( real(b') ) + 1i*cov( imag(b') ) )/100;
        end   
        
end


switch O.solver.NoiseCov_method
    
    case 1
        COV     =   eye( size(COV) ) * mean( diag(COV) );
        
    case 2
        COV     =   diag( diag( COV ) );

end


return