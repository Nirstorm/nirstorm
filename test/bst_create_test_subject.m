function [subject_name, sSubject, iSubject] = bst_create_test_subject()
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
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

subject_name = 'test_subject';
[sSubject, iSubject] = db_add_subject(subject_name);

sTemplates = bst_get('AnatomyDefaults');
iTemplate = strcmpi('Colin27_4NIRS_lowres', {sTemplates.Name});
if ~any(iTemplate)
    error('Template Colin27_4NIRS_lowres not found');
end
db_set_template(iSubject, sTemplates(iTemplate), 0);
db_save();
sSubject = bst_get('Subject', subject_name);
end
