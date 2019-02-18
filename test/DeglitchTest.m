classdef DeglitchTest < matlab.unittest.TestCase
     
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

        function test_deglitch(testCase)
            
            % Request utest data -> time-series + ground-truth marking
            deglitch_fn = utest_request_data({'artifacts','glitch_data.mat'});
            deglitch_data = load(deglitch_fn);
            
            % Import in brainstorm
            % subject_name = bst_create_test_subject('');
            deglitched_signal = process_nst_deglitch('Compute', deglitch_data.nirs_signal);
            if 0
                for ipos=1:size(deglitch_data.nirs_signal, 1)
                    figure(); hold on;
                    plot(deglitch_data.nirs_signal(ipos, :), 'b', 'LineWidth', 3);
                    plot(deglitched_signal(ipos, :), 'r');
                end
            end
            
            % Check result
            nb_samples = size(deglitch_data.nirs_signal, 2);
            for ipos=1:size(deglitch_data.nirs_signal, 1)
                gflags = deglitch_data.glitch_flags(ipos, :);
                gflags([1 end]) = 0;
                signal_clean = deglitch_data.nirs_signal(ipos, :);
                signal_clean(1) = mean(signal_clean(2:3));
                signal_clean(end) = mean(signal_clean(nb_samples-2:nb_samples-1));
                glitches_idx = find(gflags);
                signal_clean(gflags) = (signal_clean(glitches_idx-1) + signal_clean(glitches_idx+1)) / 2;
                
                if 0
                    figure(); hold on;
                    plot(deglitch_data.nirs_signal(ipos, :), 'k', 'LineWidth', 3);
                    plot(signal_clean, 'b', 'LineWidth', 2);
                    plot(deglitched_signal(ipos, :), 'r');
                end
                
                testCase.assertTrue(all_close(signal_clean, deglitched_signal(ipos, :)));
            end        
        end
        
        function test_detection(testCase)
            
            % Request utest data -> time-series + ground-truth marking
            deglitch_fn = utest_request_data({'artifacts','glitch_data.mat'});
            deglitch_data = load(deglitch_fn);
            
            glitch_flags = process_nst_deglitch('detect_glitches', deglitch_data.nirs_signal);
            
            % Check against ground-truth markings
            testCase.assertTrue(all(glitch_flags(:)==deglitch_data.glitch_flags(:)));
        end
       
        
    end
end
