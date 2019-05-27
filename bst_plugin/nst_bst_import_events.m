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



STRUCT = 1;
FILE = 2;
if isstruct(sFile)
    obj_type = STRUCT;
else
    assert(ischar(sFile));

    sFile_fn = sFile;
    obj_type = FILE;
    
    % fileType = file_gettype( fileName )
    % See process_nst_import_csv_events for proper event importation and
    % loading of sFile struct from file

    
    sFile = in_fopen_bstmat(DataFile);
end

[root, bfn, ext] = fileparts(event_fn);
switch ext
    case '.mat'
        sFile = import_events(sFile, [], event_fn, 'BST');
    otherwise
        error(sprintf('Unsupported extension for event file: "%s"', ext));
end

if obj_type == FILE
    bst_save(sFile_fn, sFile);
end

end