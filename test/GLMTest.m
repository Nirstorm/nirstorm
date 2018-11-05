classdef GLMTest < matlab.unittest.TestCase
    
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
        
        function test_simulation(testCase)
            
            % Simulate cortical signals -> dHb
            sRaw = bst_create_dummy_data_16x16(testCase.tmp_dir);
            
            sRaw = bst_create_nirs_data('simulated_evoked_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                        
            % Use head model to project signals in channel space -> dOD
            
            % Add some measurement noise
            
            % Run projection
            
            % Run surface GLM
            
        end
        
    end
    
end

function sFile = bst_create_dummy_data_16x16(tmp_dir)

%% Retrieve data
repo_url = nst_get_repository_url();
data_fns = nst_request_files({{'unittest','lesca_data','dummy_frontal_16x16.nirs'}, ...
                              {'unittest','lesca_data','optodes_frontal_16x16_Colin27_4NIRS.txt'}, ...
                              {'unittest','lesca_data','headmodel_optodes_frontal_16x16_Colin27_4NIRS.mat'}}, ...
                              1, repo_url);
nirs_fn = fullfile(tmp_dir, 'dummy_frontal_16x16.nirs');
copyfile(data_fns{1}, nirs_fn);
copyfile(data_fns{2}, fullfile(tmp_dir, 'optodes.txt'));

%% Import data in brainstorm
bst_create_test_subject();
sDummy = utest_import_nirs_in_bst(nirs_fn, 0);

%% Inject headmodel
[condition_dir, bfn, ext] = fileparts(file_fullpath(sDummy.FileName));
headmodel_fn = fullfile(condition_dir, 'headmodel_nirs_mcx_fluence.mat');
copyfile(data_fns{3},  headmodel_fn);

sStudy = bst_get('Study', sDummy.iStudy);
sStudy.iHeadModel = 1;
sStudy.HeadModel(end+1) = db_template('Headmodel');
sStudy.HeadModel(end).FileName = file_short(headmodel_fn);
sStudy.HeadModel(end).Comment = 'NIRS head model (all pairs)';
sStudy.HeadModel(end).HeadModelType = 'surface';
bst_set('Study', sDummy.iStudy, sStudy);
db_save();

%% Create scout for activating region

end