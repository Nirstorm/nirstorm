function sFileUpdated = nst_bst_import_events(sFile, events_fn)
% Wrapper around bst functions to import events from given file event_fn into 
% given sFile (file or struct), in a backward compatible way.
% Load according to format based on extension
%
%  Inputs:
%      - sFile (struct or str): data to import into. 
%              either a brainstorm data structure or a mat file containing
%              the data structure
%      - events_fn (str): event file
%              Event information. Suported format: .mat
%
%  Outputs:
%      - sFile (struct or str): same type as sFile intput
%
% TODO: utest with struct input

STRUCT = 1;
FILE = 2;
if isstruct(sFile)
    obj_type = STRUCT;
else
    assert(ischar(sFile));

    sFile_fn = sFile;
    obj_type = FILE;

    isRaw = ~isempty(strfind(sFile_fn, '_0raw'));
    if isRaw
        DataMat = in_bst_data(sFile_fn, 'F');
        sFile = DataMat.F;
    else
        DataMat = in_fopen_bstmat(sFile_fn);
        sFile = DataMat;
    end    
end

[root, bfn, ext] = fileparts(events_fn);
switch ext
    case '.mat'
        sFileUpdated = import_events(sFile, [], events_fn, 'BST');
    otherwise
        error(sprintf('Unsupported extension for event file: "%s"', ext));
end

if obj_type == FILE   
    % Report changes in .mat structure
    if isRaw
        DataMat.F = sFileUpdated;
    else
        DataMat.Events = sFileUpdated.events;
    end
    % Save file definition
    bst_save(file_fullpath(sFile_fn), DataMat, 'v6', 1);
    
    sFileUpdated = sFile_fn;
end

end
