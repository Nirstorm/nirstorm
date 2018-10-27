function file_names = nst_get_bst_func_files(subject_name, condition_name, item_name, data_types, protocol_name)

if nargin < 4
    data_types = nst_get_bst_data_fields();
else
    if ischar(data_types)
        data_types = {data_types};
    end
end

if nargin >= 5
    gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', protocol_name));
end

sSubject = bst_get('Subject', subject_name, 1);
sStudy = bst_get('StudyWithSubject', sSubject.FileName);

file_names = {};

for i_study=1:length(sStudy)
    for i_cond=1:length(sStudy(i_study).Condition)
        if strcmp(sStudy(i_study).Condition{i_cond}, condition_name)
            for i_data_type=1:length(data_types)
                data_type = data_types{i_data_type};
                for i_data=1:length(sStudy(i_study).(data_type))
                    if strcmp(sStudy(i_study).(data_type)(i_data).Comment, item_name)
                        file_names{end+1} = sStudy(i_study).(data_type)(i_data).FileName;
                    end
                end
            end
        end
    end
end

% if isempty(iData)
%     disp(['Cannot find data "' item_name ...
%            '" for subject "' acq_tag '" and condition "' condition_name '"']);
%     return 
% end
% if ~isempty(iData)
%     sFile = bst_process('GetInputStruct', sStudy(iStudy).(data_type)(iData).FileName);
% end

if length(file_names) == 1
    file_names = file_names{1};
end

end


function file_types = nst_get_bst_data_fields()

file_types = {'HeadModel', 'Channel', 'Data', 'Stat', 'Image', 'Matrix', 'Result', 'Timefreq'};

end