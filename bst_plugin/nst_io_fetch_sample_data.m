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
    otherwise
        error(['Unknown data set: ' data_label]);
end

end