function [stand_alone, process, install, bst_init] = be_check_caller()
% Check is main script is called from BST environment

stand_alone = 1;
process     = 0;
install     = 0;
bst_init    = 0;

% Check caller
caller  =   dbstack(1);
isglob  =   whos('global', 'GlobalData');

if any(is_bst_function({caller.name})) && ~isempty(isglob)
    stand_alone =  0;    
end

if any( cell2mat(strfind( {caller.name}, 'process_' )) )
    process     =  1;
end

if numel(caller)>1 && strcmp( caller(2).name, 'be_look_for_brainentropy' )
    install     =   1;
end

if any(strcmp( {caller.name}, 'bst_startup' )) 
    bst_init    =  1;    
end

return


function [status] = is_bst_function(fcns)

status  =   [];
try
    % This will fail if the call to BEst is stand alone
    bstDir  =   bst_get('BrainstormHomeDir');
    for ii  =   1 : numel(fcns)
    
        % Check current fcn
        try
            fullP   =   which(fcns{ii});
            status  =   [status ~isempty(strfind(fullP, bstDir))];
        end
    end
end

return