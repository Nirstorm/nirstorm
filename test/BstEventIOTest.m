classdef BstEventIOTest < matlab.unittest.TestCase
    
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

 
        function test_load_save_current_version(testCase)
            bst_events = utest_get_test_bst_events();
            event_fn = fullfile(testCase.tmp_dir, 'bst_events.mat');
            nst_bst_save_events(bst_events, event_fn);
            testCase.assertTrue(exist(event_fn, 'file')==2);
            loaded_events = nst_bst_load_events(event_fn);
            assert_struct_are_equal(bst_events, loaded_events);
        end
        
        function test_load_save_backward_compatibility(testCase)
            bst_events = utest_get_test_bst_events();
            % Dump event file and check that it's available online in the
            % unittest nirstorm data repository
            forge_current_utest_data_file(bst_events);
            
            % Loop through all available previous versions of dumped event
            % file: load them and check that they are still consistent with
            % current test events.
            utest_data_dir = get_utest_data_dir();
            flist = dir(utest_data_dir);
            for ifile=1:length(flist)
                if ~strcmp(flist(ifile).name, '.') && ~strcmp(flist(ifile).name, '..') 
                    loaded_events = nst_bst_load_events(fullfile(utest_data_dir, flist(ifile).name));
                    assert_struct_are_equal(bst_events, loaded_events);
                end
            end
        end
        
        % Error tests: non-available channels, incompatible time axes, 
        % Assert warning is issued when rounding of events make them
        % differ from what was stored on the drive.
    end
end

function forge_current_utest_data_file(bst_events)
% Save given bst events to file tagged with current bst version, if not
% already available.
% Insure that this file is available online.

bst_version = bst_get('Version');
event_bfn = sprintf('bst_events_v%s.mat', bst_version.Version);
event_fn = fullfile(get_utest_data_dir(), event_bfn);
if exist(event_fn, 'file') ~= 2
    nst_bst_save_events(bst_events, event_fn);
end
assert(exist(event_fn, 'file')==2);

event_fn_url = [get_utest_data_url() '/' event_bfn];
% Check if remote file exists:
jurl = java.net.URL(event_fn_url);
conn = openConnection(jurl);
conn.setConnectTimeout(15000);
status = getResponseCode(conn);
if status == 404
    error('Remote utest event data file not found: %s', event_fn_url);
end
end

function utest_data_dir = get_utest_data_dir()
sdirs = get_utest_data_subdirs();
utest_data_dir = fullfile(nst_get_local_user_dir(), sdirs{:});
if ~exist(utest_data_dir, 'dir')
    mkdir(utest_data_dir);
end
end

function utest_data_url = get_utest_data_url()
sdirs = get_utest_data_subdirs();
utest_data_url = [nst_get_repository_url() '/' strjoin(sdirs, '/')];
end

function sdirs = get_utest_data_subdirs()
sdirs = {'unittest', 'bst_events'};
end

function assert_struct_are_equal(s1, s2)

assert(length(s1) == length(s2));
f1 = fieldnames(s1);
f2 = fieldnames(s2);
assert(length(f1) == length(f2));
assert(all(ismember(f1, f2)));
for i_item=1:length(s1)
    for ifield=1:length(f1)
        f = f1{ifield};
        if ~isequaln(s1(i_item).(f), s2(i_item).(f))
            error('Unequal values for field %s', f);
        end
    end
end
end