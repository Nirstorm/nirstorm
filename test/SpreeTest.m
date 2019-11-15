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

            spree_result_chan = bst_process('CallProcess', 'process_nst_spree', nirs_chan_fn, [], ...
                                       'lf_cutoff', 0.01, ...
                                       'trim_start', 0, ...
                                       'stim_events', strjoin({activity.condition_name}, ','));
            check_chan_activity_detection(spree_result_chan.activ_chans, activity);
            
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

% Request activity test data

% Import data in bst

% Retrieve file names

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