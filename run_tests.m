function run_tests()
%TODO: check for matlab version > R2013b
%TODO: check that nirstorm is actually installed

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite
import matlab.unittest.plugins.StopOnFailuresPlugin
import matlab.unittest.plugins.CodeCoveragePlugin

%% Sort out tests in different test suites:
%   - tests for source tools (eg installation)
%   - tests for installed package (eg brainstorm processes)
test_scripts = dir(fullfile('test', '*.m'));
iss = 1;
ips = 1;
for iscript=1:length(test_scripts)
    test_fn = fullfile('test', test_scripts(iscript).name);
    tests = TestSuite.fromFile(test_fn);
    if ~isempty(strfind(test_fn, 'SourceTest'))
        source_suite(iss:(iss+length(tests)-1)) = tests;
        iss = iss + length(tests);
    else
        package_suite(ips:(ips+length(tests)-1)) = tests;
        ips = ips + length(tests);
    end
end

%% Configure & run test runner for source tools
runner = TestRunner.withTextOutput;
runner.addPlugin(StopOnFailuresPlugin('IncludingAssumptionFailures', true));
runner.addPlugin(CodeCoveragePlugin.forFolder('dist_tools'));
result = runner.run(source_suite);

%% Configure & run test runner for installed package tools
runner = TestRunner.withTextOutput;
runner.addPlugin(StopOnFailuresPlugin('IncludingAssumptionFailures', true));
runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(bst_get('BrainstormUserDir'), 'process')));
% result = runner.run(package_suite);
end