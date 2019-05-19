classdef SurfTplPipelineImportTest < matlab.unittest.TestCase
    
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
        
        function test_import_subjects(testCase)
            % TODO
%             nirs_fns = simulate_nirs();
%             options = nst_ppl_surface_template_V1('setup');
%             options.moco.export_dir = fullfile(testCase.tmp_dir, 'moco');
%             sFilesRaw = nst_ppl_surface_template_V1('import', {'data1.nirs', 'data2.nirs'}, {'subj1', 'subj2'});


        end

        function test_import_events(testCase)

        end

        function test_import_bad_channels(testCase)

        end
        
    end
end
