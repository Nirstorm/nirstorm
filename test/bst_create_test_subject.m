function [subject_name, sSubject, iSubject] = bst_create_test_subject(template_name, use_default_anatomy)

if nargin < 1
    template_name = 'Colin27_4NIRS_lowres';
end

if nargin < 2
    use_default_anatomy = 0;
end

bst_create_test_protocol(use_default_anatomy);

% Add test subject
subject_name = 'test_subject';
[sSubject, iSubject] = db_add_subject(subject_name);

if ~isempty(template_name)
    % Set anatomy to template Colin27 4NIRS - low res version
    sTemplates = bst_get('AnatomyDefaults');
    iTemplate = strcmpi(template_name, {sTemplates.Name});
    if ~any(iTemplate)
        error('Template %s not found', template_name);
    end
    db_set_template(iSubject, sTemplates(iTemplate), 0);
    db_save();
end
% Set default cortical surface to low resolution mid
sSubject = bst_get('Subject', subject_name);
db_surface_default( iSubject, 'Cortex', find(strcmp({sSubject.Surface.Comment}, 'mid_lowres')));
panel_protocols('RepaintTree');
end
