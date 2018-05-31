classdef EVTImportTest < matlab.unittest.TestCase
    
    properties
        tmp_dir
    end
    
    
    % This script aimed to test two functionality :
    % 
    %
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
        function test_evt_event_import_events(testCase)
            global GlobalData;
            
            sFiles = {...
            'Subject01/NIRSCOMP_S01-2017-07-12_001_Hb/data_hb_180423_1149_bl_detrend_band.mat'};
            RawFiles = {...
            '/lena13/home_users/users/ext-e.delaire/Desktop/test_event/test_event/NIRSCOMP_S01-2017-07-12_001.evt'};
            
            sFiles_processed = bst_process('CallProcess', ...
                'process_nst_import_evt_events', ...
                sFiles, [],  ...
                'evtfile',{RawFiles{1}, ''}, ...
                'label_cols','ConditionA,ConditionB,ConditionC', ...
                'confirm_importation', 0);
                           
                           
            assert(~isempty(sFiles_processed));
            events = get_events(sFiles_processed);
            
            assert(length(events) == 3);
            assert(size(events(1).samples, 1) == 2);
            
            %assert(abs(events(strcmp({events.label}, 'DMNirs')).times(1, 1) - 1482333540.553) <= nirs_dt);
            %assert(abs(events(strcmp({events.label}, 'DMNirs')).times(2, 1) - 1482333570.823) <= nirs_dt);
            utest_reset_bst();      
        end
        
        
        function test_evt_event_merge(testCase)
            global GlobalData;
            
            sFiles = {...
            'Subject01/NIRSCOMP_S01-2017-07-12_001_Hb/data_hb_180423_1149_bl_detrend_band.mat'};
            RawFiles = {...
            '/lena13/home_users/users/ext-e.delaire/Desktop/test_event/test_event/NIRSCOMP_S01-2017-07-12_001.evt'};
            
            sFiles_processed = bst_process('CallProcess', ...
                'process_nst_import_evt_events', ...
                sFiles, [],  ...
                'evtfile',{RawFiles{1}, ''}, ...
                'label_cols','ConditionA,ConditionB,ConditionB', ...
                'confirm_importation', 0);
                           
                           
            assert(~isempty(sFiles_processed));
            events = get_events(sFiles_processed);
            
            assert(length(events) == 2);
            assert(size(events(1).samples, 1) == 2);
            utest_reset_bst();      
        end
        
        
    end
end

function events = get_events(sFile)
DataMat = in_bst_data(sFile.FileName, 'F');
events = DataMat.F.events;
end
