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
            rmdir(fullfile(nst_get_local_user_dir(),'unittest',testCase.tmp_dir), 's');
            utest_clean_bst();
        end
    end
    
    methods(Test)
        function test_evt_event_import_events(testCase)
            global GlobalData;
            
            utest_dir=fullfile(nst_get_local_user_dir(),'unittest');
            if ~exist(fullfile(utest_dir,testCase.tmp_dir), 'dir')
                mkdir(fullfile(utest_dir,testCase.tmp_dir));
            end
            
            events=generate_fake_event(10,2,150);
            
            dlmwrite(fullfile(utest_dir,testCase.tmp_dir,'event.evt'),events,'delimiter','\t')
            %bst_save(fullfile(utest_dir,testCase.tmp_dir,'event.evt'),events);
            %bst_save(fullfile(utest_dir,'event.evt') ,events);

            
            nirs = load(fullfile(utest_dir,'dummy.nirs'), '-mat');
            nirs_dt = diff(nirs.t(1:2));
            
            
            sFiles_dummy = utest_import_nirs_in_bst(fullfile(utest_dir,'dummy.nirs'));

            sFiles_processed = bst_process('CallProcess', ...
                'process_nst_import_evt_events', ...
                sFiles_dummy, [],  ...
                'evtfile',{ fullfile(utest_dir,testCase.tmp_dir,'event.evt')}, ...
                'label_cols','ConditionA,ConditionB', ...
                'last_event',0 ,...
                'confirm_importation', 0);
                           
                           
            assert(~isempty(sFiles_processed));
            events = get_events(sFiles_processed);
            
            assert(length(events) == 2);
            assert(size(events(1).samples, 1) == 2);
            
            
            assert( strcmp( events(1).label,'ConditionB'));
            assert(  size(events(1).epochs,2) == 9 ); 
            % The last event is not here because we don't know his duration
            % at it is the last of the file

            assert( strcmp( events(2).label,'ConditionA'));
            assert(  size(events(2).epochs,2) == 10 );

            
            assert(abs( events(strcmp({events.label}, 'ConditionA')).times(1, 1)) - 6.04000000000000 <= nirs_dt);
            assert(abs( events(strcmp({events.label}, 'ConditionA')).times(2, 1)) - 6.04000000000000 <= nirs_dt);
            
            utest_reset_bst();      
        end
        
        
        function test_evt_event_merge(testCase)
            global GlobalData;
            
                        utest_dir=fullfile(nst_get_local_user_dir(),'unittest');
            if ~exist(fullfile(utest_dir,testCase.tmp_dir), 'dir')
                mkdir(fullfile(utest_dir,testCase.tmp_dir));
            end
            
            events=generate_fake_event(10,4,150);
            
            dlmwrite(fullfile(utest_dir,testCase.tmp_dir,'event.evt'),events,'delimiter','\t')
            %bst_save(fullfile(utest_dir,testCase.tmp_dir,'event.evt'),events);
            %bst_save(fullfile(utest_dir,'event.evt') ,events);

            
            nirs = load(fullfile(utest_dir,'dummy.nirs'), '-mat');
            nirs_dt = diff(nirs.t(1:2));
            
            
            sFiles_dummy = utest_import_nirs_in_bst(fullfile(utest_dir,'dummy.nirs'));

            sFiles_processed = bst_process('CallProcess', ...
                'process_nst_import_evt_events', ...
                sFiles_dummy, [],  ...
                'evtfile',{ fullfile(utest_dir,testCase.tmp_dir,'event.evt')}, ...
                'label_cols','ConditionA,ConditionA,ConditionB,ConditionC', ...
                'last_event',0 ,...
                'confirm_importation', 0);
            
            assert(~isempty(sFiles_processed));
            events = get_events(sFiles_processed);
            
            assert(length(events) == 3);
            assert(events(strcmp({events.label}, 'ConditionA')).times(2, 1) - events(strcmp({events.label}, 'ConditionB')).times(1, 1) - 12 <= nirs_dt);
            assert(events(strcmp({events.label}, 'ConditionB')).times(2, 1) - events(strcmp({events.label}, 'ConditionB')).times(1, 1) - 6 <= nirs_dt);
            assert(events(strcmp({events.label}, 'ConditionC')).times(2, 1) - events(strcmp({events.label}, 'ConditionC')).times(1, 1) - 6 <= nirs_dt);

            
            utest_reset_bst();      
        end
        
        
    end
end

function events = get_events(sFile)
DataMat = in_bst_data(sFile.FileName, 'F');
events = DataMat.F.events;
end


% Generate a file with n_unique events which repete n times. 
% Each event last duration sample. 

function fake_event=generate_fake_event(n,unique_event,duration)

fake_event=zeros(n*unique_event,9);
sample=1;

for i=1:(n*unique_event) 
        fake_event(i,1)=sample;
        event_code=dec2bin(1+mod(i,unique_event),8);
        for j=2:9
            fake_event(i,j)=str2num(event_code(j-1));          
        end
        sample=sample+duration;
end
end