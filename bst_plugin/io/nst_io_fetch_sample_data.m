function varargout = nst_io_fetch_sample_data(data_label, confirm_download)
% Return file names of sample data. Download them if necessary.
%
% [NIRS_FILES, SUBJECT_NAMES] = NST_IO_FETCH_SAMPLE_DATA('template_group_tapping', CONFIRM_DOWNLOAD)
%

if nargin < 2
    confirm_download = 1;
end

switch data_label
    case 'template_group_tapping'
        subject_ids = 1:10;
        subject_names = arrayfun(@(ii) sprintf('Subject%02d', ii), subject_ids, 'UniformOutput', 0);
        nb_subjects = length(subject_ids);

        data_to_fetch = cell(1, nb_subjects * 2);
        % Make list of nirs files to fetch:
        % sample_data/template_group_tapping/Subject*/S*_tapping.nirs
        for subject_id=subject_ids
            data_to_fetch{subject_id} = {'sample_data', 'template_group_tapping', ...
                                         subject_names{subject_id}, ... 
                                         sprintf('S%02d_tapping.nirs', subject_id)};
        end
        % Make list of optode files to fetch:
        % sample_data/template_group_tapping/Subject*/optodes.txt
        for subject_id=subject_ids
            data_to_fetch{nb_subjects + subject_id} = {'sample_data', 'template_group_tapping', ...
                                                       subject_names{subject_id}, 'optodes.txt'};
        end
        data_fns = nst_request_files(data_to_fetch, confirm_download, nst_get_repository_url(), 1e6);
        varargout{1} = data_fns(1:nb_subjects);
        varargout{2} = subject_names;
        
    case 'group_tapping'
        subject_names = {'S01' ,'S02','S03',  ...
                         'S04', 'S05',        ...
                         'S07', 'S08', 'S09', ...
                         'S10', 'S11'};
             
        nb_subjects = length(subject_names);
        data_to_fetch{3*nb_subjects}={''};

        for i=1:nb_subjects
            data_to_fetch{i}={'Tapping',subject_names{i}, 'data', ['subject_' subject_names{i}(2:3) '.nirs'] };
            data_to_fetch{i+nb_subjects}={'Tapping',subject_names{i}, 'data', 'optodes.txt' };
            data_to_fetch{i+2*nb_subjects}={'Tapping',subject_names{i}, 'data', 'headpoints' };  
        end
        
        data_fns = nst_request_files(data_to_fetch, confirm_download, nst_get_repository_url(), 1e6);                    
        varargout{1} = data_fns(:);
        varargout{2} = subject_names;                
    case 'group_tapping_with_anatomy'
        subject_names = {'S01' ,'S02','S03',  ...
                         'S04', 'S05',        ...
                         'S07', 'S08', 'S09', ...
                         'S10', 'S11'};
             
        nb_subjects = length(subject_names);
        data_to_fetch{4*nb_subjects}={''};

        for i=1:nb_subjects
            data_to_fetch{i}={'Tapping',subject_names{i}, 'anatomy' };
            data_to_fetch{i+nb_subjects}={'Tapping',subject_names{i}, 'data', ['subject_' subject_names{i}(2:3) '.nirs'] };
            data_to_fetch{i+2*nb_subjects}={'Tapping',subject_names{i}, 'data', 'optodes.txt' };
            data_to_fetch{i+3*nb_subjects}={'Tapping',subject_names{i}, 'data', 'headpoints' };
        end
        
        data_fns = nst_request_files(data_to_fetch, confirm_download, nst_get_repository_url(), 1e6);                    
        varargout{1} = data_fns(:);
        varargout{2} = subject_names; 
        
   
    otherwise
        error(['Unknown data set: ' data_label]);
end

end