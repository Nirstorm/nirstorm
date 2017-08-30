function run_tests(to_run, stop_on_error)
% Run unit test suites from files located in './test'
% Test scenarios to run are specificied by given to_run:
%    - package: tests related to functions in bst_plugin
%    - source: tests of tools used on package sources (eg dist_tools)
%    - scripts: execute tutorial script (download data if not available).
%               WARNING: script tests take time!!
% 
install_error_msg = sprintf(['To run unit tests, nirstorm debug functions must be installed.\n'...
                            ' Use nst_install(''copy'', ''debug'') or ' ...
                            'nst_install(''copy'', ''debug'') (linux only).\n'...
                            'WARNING: these functions override brainstorm behavior.']);
try
    if ~exist(fullfile(bst_get('BrainstormUserDir'), 'process', 'bst_error'), 'file');
        error(install_error_msg);
    end
catch
xà    error(install_error_msg);
end
%TODO: check for matlab version > R2013b
%TODO: check that nirstorm is actually installed

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.plugins.CodeCoveragePlugin

if nargin < 1
    to_run = {'package', 'source'};
else
    if ischar(to_run) 
        if strcmp(to_run, 'all')
            to_run = {'package', 'source', 'scripts'};
        else
            to_run = {'package'};
            end
    end
        
    assert(all(ismember(to_run, {'package', 'source', 'scripts'})));
end

if nargin < 2
    stop_on_error = 0;
end

%% Sort out tests in different test suites:
%   - tests for source tools (eg installation)
%   - tests for installed package (eg brainstorm processes)
test_scripts = dir(fullfile('test', '*.m'));
iss = 1;
ips = 1;
ics = 1;
for iscript=1:length(test_scripts)
    test_fn = fullfile('test', test_scripts(iscript).name);
    if ~isempty(strfind(test_fn, 'Test.m'))
        tests = TestSuite.fromFile(test_fn);
        if ~isempty(strfind(test_fn, 'SourceTest'))
            source_suite(iss:(iss+length(tests)-1)) = tests;
            iss = iss + length(tests);
        elseif ~isempty(strfind(test_fn, 'ScriptTest'))
            script_suite(ics:(ics+length(tests)-1)) = tests;
            ics = ics + length(tests);
        else
            package_suite(ips:(ips+length(tests)-1)) = tests;
            ips = ips + length(tests);
        end
    end
end

if ismember('package', to_run) && ~isempty(package_suite)
    %% Configure & run test runner for installed package tools
    runner = TestRunner.withTextOutput;
    if stop_on_error
        runner.addPlugin(StopOnFailuresPlugin('IncludingAssumptionFailures', true));
    end
    runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(bst_get('BrainstormUserDir'), 'process')));
    result = runner.run(package_suite);
end

if ismember('source', to_run)
    %% Configure & run test runner for source tools
    runner = TestRunner.withTextOutput;
    if stop_on_error
        runner.addPlugin(StopOnFailuresPlugin('IncludingAssumptionFailures', true));
    end
    addpath('dist_tools');
    runner.addPlugin(CodeCoveragePlugin.forFolder('dist_tools'));
    result = runner.run(source_suite);
end

if ismember('scripts', to_run)
    %% Configure & run test runner for source tools
    runner = TestRunner.withTextOutput;
    if stop_on_error
        runner.addPlugin(StopOnFailuresPlugin('IncludingAssumptionFailures', true));
    end
    result = runner.run(script_suite);
end

end