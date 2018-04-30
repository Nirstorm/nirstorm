function [a,b,c] = be_fileparts(fileName)
% This function returns the root folder=a, file name=b and extension=c of 
% the file named FILENAME 

% Prepare output
a   =   fileName;
b   =   '';
c   =   '';

% Get info
idPnt   =   strfind( fileName, '.' );
if ~isempty(idPnt)
    c   =   fileName(idPnt(end):end);
    fileName(idPnt(end):end)    =   [];
end

idSep   =   strfind( fileName, filesep );
if ~isempty(idSep) 
    a   =   fileName(1:idSep(end)-1);
    b   =   fileName(idSep(end)+1 : end);
end
    
return

