classdef SpreeTest < matlab.unittest.TestCase
    
    properties
        tmp_dir
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            tmpd = tempname;
            mkdir(tmpd);
            testCase.tmp_dir = tmpd;
            utest_bst_setup();
        end
    end
    
    methods(TestMethodTeardown)
        function tear_down(testCase)
            rmdir(testCase.tmp_dir, 's');
            utest_clean_bst();
        end
    end
    
    methods(Test)
        function test_detection(testCase)
            % Load forged data to test activity detection
            % activity is a struct array with fields:
            %  - condition_name (str)
            %  - activ_chans: cell array of channel indexes (int) that are
            %                 activating for this condition
            %  - activ_scout: bst scout gathering all cortical vertices
            %                 that are activating for this condition
            [nirs_chan_fn, nirs_cortex_hbo, nirs_cortex_hbr, activity] = load_activity_test_subject();

            
            if 1
                % Requires figtk
                fig_options.DefaultFigureVisible = 'on';
                fig_options.fig_height = 5;
                fig_options.fig_width = 6;
                fig_options.DefaultAxesFontSize = 22;
                fig_options.DefaultTextFontSize = 22;
                figtk_setup(fig_options);
                
                output_fig_dir = '/tmp/nst_spree_utest/';
            else
                output_fig_dir = '';
            end
            
            spree_result_chan = bst_process('CallProcess', 'process_nst_spree', nirs_chan_fn, [], ...
                                            'stim_events', 'stim2', ...
                                            'nb_iterations', 200, ...
                                            'save_fit', 1, ...,
                                            'output_fig_dir', output_fig_dir, ...
                                            'save_full_fitted_model', 1);
                         
            % Only check S1D1
            check_response_estimation(spree_result_chan);                                
            check_chan_activity_detection(spree_result_chan.activ_chans, activity);
            
            % Only check S2D2
            
            
            spree_result_cortex_hbo = bst_process('CallProcess', 'process_nst_spree', ...
                                                  nirs_cortex_hbo, [], ...
                                                  'lf_cutoff', 0.01, ...
                                                  'trim_start', 0, ...
                                                  'stim_events', strjoin({activity.condition_name}, ','));
            check_cortex_activity_detection(spree_result_cortex_hbo.activ_scouts, activity);
        end
    end
end



function [chan_fn, cortex_hbo, cortex_hbr, activity] = load_activity_test_subject()

% Setup utest protocol
use_default_anatomy = 0;
bst_create_test_protocol(use_default_anatomy);

% Request activity test data
repo_url = nst_get_repository_url();
%TODO: add ground-truth HRF in forged subject data
data_fns = nst_request_files({{'unittest','activity_test_subject','activity_test_subject.zip'}, ...
                              {'unittest','activity_test_subject','activity_info.mat'}}, ...
                              1, repo_url);
                          
% Import data in bst
import_subject(data_fns{1});

% Retrieve file names
% TODO: rename test_subject to activity_test_subject
%       and condition to "evoked_activity_high_SNR"
chan_fn = nst_get_bst_func_files('test_subject', 'test', 'Hb');
cortex_hbo = nst_get_bst_func_files('test_subject', 'test', 'HbO cortex');
cortex_hbr = nst_get_bst_func_files('test_subject', 'test', 'HbR cortex');

activity_info_data = load(data_fns{2}, '-mat');
activity = activity_info_data.activity_info;
end

function check_chan_activity_detection()

% Sensitivity check:
% All activ channels are actually detected


% Specificity check:
% All inactiv channels are not detected

end

function check_cortex_activity_detection()

% Sensitivity check:
% All activ vertices are actually detected


% Specificity check:
% All inactiv vertices are not detected

end
