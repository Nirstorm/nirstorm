function signal= Constant(t,value) 
% Constant : return a constant function.
    
    %signal=zeros( size(t_vect) );  
    if nargin <  2, value = 1; end

    assert( isvector(t)  )
    signal = value * ones(numel(t),1);
    
end

