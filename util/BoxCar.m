function signal= BoxCar(t_vect,lag,duration)
% BoxCar : apply the box car function over t_vect
    
    signal=zeros( size(t_vect) );
    i=1;
    
    for t=t_vect
        if( t >= lag && t < lag + duration ) 
            signal(i)=1;
        else
            signal(i)=0;
        end;
        i=i+1;
    end
end

