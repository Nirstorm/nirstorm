function sFileUpdated = nst_bst_import_channel_flags(sFile,  channel_flags_fn)
%NST_BST_EXPORT_CHANNEL_FLAGS Save channel flags from given sFile 
%                             into given channel_flags_fn.
% TODO: also save list of channels to always insure that chan flags array
% is aligned with channel info
%
%  Inputs:
%      - sFile (struct or str): data to import into. 
%              either a brainstorm data structure or a mat file containing
%              the data structure
%      - events_fn (str): channel flags file
%              Channel flag information. Suported format: .mat
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

[root, bfn, ext] = fileparts(channel_flags_fn);
sFileUpdated = sFile;
switch ext
    case '.mat'
        chan_flags_data = load(channel_flags_fn);
        sFileUpdated.ChannelFlag = chan_flags_data.channel_flags;
        %TODO: save channel flags properly -> should be in channel header
        %file, not data file!
    otherwise
        error(sprintf('Unsupported extension for channel file: "%s"', ext));
end

if obj_type == FILE   
    % Report changes in .mat structure
    DataMat.ChannelFlag = sFileUpdated.ChannelFlag;
    if isRaw
        DataMat.F.channelflag = chan_flags_data.channel_flags;
    end
    % Save file definition
    bst_save(file_fullpath(sFile_fn), DataMat, 'v6', 1);
    
    sFileUpdated = sFile_fn;
end

end
