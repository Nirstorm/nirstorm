classdef SurfTplPipelineTest < matlab.unittest.TestCase
    
    properties
        tmp_dir
        ppl_tag
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            tmpd = tempname;
            mkdir(tmpd);
            testCase.tmp_dir = tmpd;
            utest_bst_setup();
            
            ProtocolName = 'nst_utest';
            %% Ensure that nst_utest protocol exists
            if isempty(bst_get('Protocol', ProtocolName))
                % Create new protocol
                gui_brainstorm('CreateProtocol', ProtocolName, 1, 0); %UseDefaultAnat=1, UseDefaultChannel=0.
            end
            
            testCase.ppl_tag = '__nspst_V1';
        end
    end
    
    methods(TestMethodTeardown)
        function tear_down(testCase)
            rmdir(testCase.tmp_dir, 's');
            utest_clean_bst();
        end
    end
    
    methods(Test)
        
        % TODO: test figure output & rewriting
        % TODO: check coverage

%         function test_import_subjects_from_scratch(testCase)
%             [nirs_fns, subject_names] = load_test_group_data();
%             sFiles = nst_ppl_surface_template_V1('import', [], nirs_fns, subject_names);
%             
%             % Check that all origin/Raw items are there
%             for isubject=1:length(subject_names)
%                 sFileRaw = nst_get_bst_func_files(subject_names{isubject}, ['origin' testCase.ppl_tag], 'Raw');
%                 testCase.assertNotEmpty(sFileRaw);
%                 testCase.assertMatches(sFiles{isubject}, sFileRaw);
%                 % Check that event group for movement artefacts is there
%                 data_event = load(file_fullpath(sFileRaw), 'Events');
%                 events = data_event.Events;
%                 testCase.assertTrue(ismember('movement_artefacts', {events.label}));
%             end
%         end

%         function test_reimportation(testCase)
%             % Do a first importation of all data, then remove one item and
%             % check that only this one is reimported while running the
%             % importation again
%             [nirs_fns, subject_names] = load_test_group_data();
%             sFiles = nst_ppl_surface_template_V1('import', [], nirs_fns, subject_names);
%             bst_process('CallProcess', 'process_delete', sFiles{2}, [], ...
%                        'target', 1);
%             [sFiles, reimported] = nst_ppl_surface_template_V1('import', [], nirs_fns, subject_names);
%             testCase.assertTrue(all(reimported==[0 1 0]));
%             sFileRaw = nst_get_bst_func_files(subject_names{2}, ['origin' testCase.ppl_tag], 'Raw');
%             testCase.assertNotEmpty(sFileRaw);
%             testCase.assertMatches(sFiles{2}, sFileRaw);
%         end


% 
        function test_save_markings(testCase)
            % Do a first importation of all data, then add some events &
            % tag some channels.
            % Run importation again and check that user-defined markings are 
            % saved.
            % TODO: add action 'save_markings' to explicitely save all
            %       events and bad channels tagging.
            
            % TODO: reimplement 
            %    - utest them with MWUC from ppl function & insure backward compatibility
            %    - see if we still need to keep nst_bst_save_events / nst_bst_load_events
            [nirs_fns, subject_names] = load_test_group_data();
            options = nst_ppl_surface_template_V1('get_options');
            event_dir = fullfile(testCase.tmp_dir, 'export_events');
            options.export_dir_events = event_dir;
            chan_flags_dir = fullfile(testCase.tmp_dir, 'export_channel_flags');
            options.export_dir_channel_flags = chan_flags_dir;
            
            sFiles = nst_ppl_surface_template_V1('import', options, nirs_fns, subject_names);
            data = load(file_fullpath(sFiles{2}));
            events = data.Events;
            i_evt_mvt = find(strcmp('movement_artefacts', {events.label}));
            testCase.assertEqual(length(i_evt_mvt), 1);
            
            % Add one movement artefacts
            events(i_evt_mvt).times = [1.5;...
                                       1.87];
            events(i_evt_mvt).epochs = [1];
            events(i_evt_mvt).reactTimes = [];
            events(i_evt_mvt).channels = {{}};
            events(i_evt_mvt).notes = {'participant sneezed'};
            
            % Append test events
            events = [events utest_get_test_bst_events()];
            
            % Tag some channels as bad:
            data.ChannelFlag(1) = -1;
            data.ChannelFlag(end) = -1;
            
            % Update bst data 
            data.Events = events;
            save(file_fullpath(sFiles{2}), '-struct', 'data'); %TODO: use bst import function instead of hacking data file
            
            % Save markings
            sFiles_mkgs = nst_ppl_surface_template_V1('save_markings', options, subject_names);
            testCase.assertMatches(sFiles_mkgs{1}, sFiles{1});
            testCase.assertMatches(sFiles_mkgs{2}, sFiles{2});
            testCase.assertMatches(sFiles_mkgs{3}, sFiles{3});
            testCase.assertEqual(length(sFiles_mkgs), length(sFiles));
            
            % Check event markings
            get_evt_fn = @(ii) fullfile(event_dir, [subject_names{ii} '_events.mat']);
            testCase.assertTrue(exist(get_evt_fn(1), 'file')==2);
            loaded_events = load(get_evt_fn(1));
            loaded_events = loaded_events.events;
            testCase.assertTrue(length(loaded_events) == 1);
            testCase.assertMatches(loaded_events(1).label, 'movement_artefacts');
            
            event_fn = get_evt_fn(2);
            testCase.assertTrue(exist(event_fn, 'file')==2);
            loaded_events = load(event_fn);
            loaded_events = loaded_events.events;
            testCase.assertTrue(isequaln(loaded_events, events));
            
            % Check bad channel markings
            get_ctag_fn = @(ii) fullfile(chan_flags_dir, [subject_names{ii} '_channel_flags.mat']);
            testCase.assertTrue(exist(get_ctag_fn(1), 'file')==2);
            
            ctag_fn = get_ctag_fn(2);
            testCase.assertTrue(exist(ctag_fn, 'file')==2);
            loaded_flags = load(ctag_fn);
            loaded_flags = loaded_flags.channel_flags;
            testCase.assertEqual(loaded_flags(1), -1);
            testCase.assertEqual(loaded_flags(end), -1);
            testCase.assertTrue(all(loaded_flags(2:end-1) == 1));
        end
% 
%         function test_reimportation_with_markings(testCase)
%             % Do a first importation of all data, then add some event & tag chans 
%             % Save event markings. Then remove item. Reimport and check 
%             % that user-defined event markings were restored

%         end
% 
%         function test_reimportation_with_markings_events_only(testCase)
%             % Do not enable saving of bad chans tags.
%             % Do a first importation of all data, then add some event & chan tags
%             % Save markings. Then remove item. Reimport and check 
%             % that user-defined event markings were restored but not chan tags.
%         end
%         
%         function test_reimportation_with_markings_bad_chans_only(testCase)
%             % Do not enable saving of events.
%             % Do a first importation of all data, then add some event & chan tags
%             % Save markings. Then remove item. Reimport and check 
%             % that chan tags were restored but not event markings.
%         end
% 
%         function test_default_pipeline(testCase)
%             % Import data set. Run pipeline with minimal options.
%             % Insure that the following outputs are produced:
%             %   - SCI
%             %   - resampled data
%             %   - dOD
%             %   - Projected data
%             %   - GLM 1st level with contrast effect maps
%             %   - GLM 2nd level with contrast t-maps
%         end
        
    end
end

function [nirs_fns, subject_names] = load_test_group_data()
[nirs_fns, subject_names] = nst_io_fetch_sample_data('template_group_tapping');
nirs_fns = nirs_fns(1:3);
subject_names = subject_names(1:3);
end