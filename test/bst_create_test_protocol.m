function bst_create_test_protocol(use_default_anatomy)

if nargin < 1
    use_default_anatomy = 0;
end

%% Ensure that nst_utest protocol exists
ProtocolName = 'nst_utest';
% Delete existing protocol
db_dir = bst_get('BrainstormDbDir');
gui_brainstorm('DeleteProtocol', ProtocolName);

nst_protocol_dir = fullfile(db_dir, ProtocolName);
if exist(nst_protocol_dir, 'dir')
    rmdir(nst_protocol_dir, 's');
end

% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, use_default_anatomy, 0);

end