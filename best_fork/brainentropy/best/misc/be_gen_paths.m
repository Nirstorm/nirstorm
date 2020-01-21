function be_gen_paths

Cpath   =   which(mfilename);

idSep   =   strfind( Cpath, filesep );
Cpath(idSep(end):end) = [];

warning('OFF')
addpath( genpath(Cpath) );        
                
return