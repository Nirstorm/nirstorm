function nst_bst_set_template_anatomy(template_name, iSubject, confirm_download)
%NST_BST_SET_TEMPLATE_ANATOMY Use given template to set the anatomy.
%   Set the anatomy of the given subject (default is 0) to the given 
%   Nirstorm template, for the currently loaded brainstorm protocol.
%
%   WARNING: It does not handle Brainstorm templates. To do that, use:
%            sTemplates = bst_get('AnatomyDefaults');
%            iTemplate = find(strcmpi(template_name, {sTemplates.Name}));
%            db_set_template(iSubject, sTemplates(iTemplate), 0);
%
%   TEMPLATE_NAME (str): name of the nirstorm template
%            See NST_CORE_GET_AVAILABLE_TEMPLATES for the description of
%            the available templates and their string identifiers
%  [ISUBJECT (int):] subject ID
%            Default is 0 (default subject)
%  [CONFIRM_DOWNLOAD (bool):] 
%            Flag to ask the user to confirm download if the template is not 
%            available locally. Default is 1.
%

if nargin < 2
    iSubject = 0;
end

if nargin < 3
    confirm_download = 1;
end

sTemplates = bst_get('AnatomyDefaults');
iTemplate = strcmpi(template_name, {sTemplates.Name});
if ~any(iTemplate)
    template_bfn = [template_name '.zip'];
    template_tmp_fn = nst_request_files({{'template', template_bfn}}, confirm_download, ...
                                        nst_get_repository_url(), 18e6);
    % Copy to .brainstorm/defaults/anatomy
    copyfile(template_tmp_fn{1}, fullfile(bst_get('BrainstormUserDir'), 'defaults', 'anatomy'));
    % Remove temporary file
    delete(template_tmp_fn{1});
end
sTemplates = bst_get('AnatomyDefaults');
iTemplate = strcmpi(template_name, {sTemplates.Name});
db_set_template(iSubject, sTemplates(iTemplate), 0);
db_save();
end

