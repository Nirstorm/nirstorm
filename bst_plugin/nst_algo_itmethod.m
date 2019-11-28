function [iterated_model]  = itmethod(model, update_variables, max_iterations, varargin) 
% itmethod - iterative method implementation that sequentially updates a
% set of variables according to the given update rule and manage history 
% tracking of variables and observables.
% Obervables are quantities computed from the updated variables via the
% given compute_observables function handle (eg compute online average).
% 
% This function can be used to implement a gradient descent, EM algorithm
% or MCMC sample scheme ...
%
% Synopsis
%   [iterated_model] = itmethod(model, update_variables, max_iterations,
%                               [var_history_pace, compute_observables,],
%                               [obs_cpt_start, obs_history_pace,]
%                               [stop_criterion,])
%
% Description
%   Starting from the initial given variables, sequentially update them
%   until max_iterations is reached or stop_criterion is met. 
%   If the function handle *compute_observable* is
%   provided then observables are also computed. History tracking is
%   defined by *var_history_pace* and *obs_history_pace* to manage the 
%   recording frequency.
%
%   Internally, the model input variable gathers all information relative 
%   to the iterative method. Ot is passed to function handles 
%   *update_variables* and *compute_observables*.
%   The model variable has the following fields:
%      variables     (struct of arrays)  % current set variables 
%                                          (defined by input)
%      iteration     (integer)           % current iteration
%      obs_iteration (integer)           % current observable iteration
%      observables   (struct of arrays)  % current set of observables
%      var_history   (struct of arrays)  % history of variables
%      obs_history   (struct of arrays)  % history of observables
%      history_iterations (map)          % keep track of recorded
%                                          iterations
%   The final state of this world variable is returned by the function
%
% Inputs ([]s are optional)
%   (struct) model
%       main data container with the following fields:
%           - (struct of arrays) variables
%              Model variables to be updated
%           - [(struct of binary flags) variables_to_estimate]
%              Indicate which variables to estimate. If a variable
%              is not estimated, its value is unchanged.
%              Should be used by given function *update_variables*
%           - [(struct) constants] 
%             Any non-changing useful data. This field is not actually
%             used here.
%           - [(struct) tmp]
%             Any preallocated quantities. This field is not actually 
%             used here.
%   (function handler) update_variables
%       function with synopsis: [model_upd] = update_variables(model)
%       where *model_upd.variables* is an updated version of the input
%       *model.variables*
%   (integer) max_iterations
%       maximum number of iterations
%   [(int) var_hist_pace]
%       Regularity of history recording for variables. 
%       If -1 then do not record (default).
%   [(function handler) compute_observables]
%       function with synopsis: [new_obs] = compute_observables(world)
%       new_obs is a struct whose fields are observable quantities (arrays)
%   [(int) obs_cpt_start]
%       Iteration index when to start computing (and recording)
%       observables.
%   [(int) obs_hist_pace]
%       Regularity of history recording for obserbables. 
%       If -1 then do not record (default).
%   [(function handler) stop_criterion]
%       if integer then it corresponds to the maximum number of iterations
%       if function handler then it should be a predicate taking the world
%           struct as argument and returning true if the iterative method
%           must terminate.
%   [(cell array of str) not_to_track]
%       List of variable or observable names for which to disable history tracking
%
% Outputs ([]s are optional)
%   (struct) iterated_model
%     The final state of the model
%     The model variable has the following fields:
%      variables     (struct of arrays)  % final values of variables 
%      iteration     (integer)           % final iteration
%      obs_iteration (integer)           % final observable iteration
%      observables   (struct of arrays)  % final observable values
%      var_history   (struct of arrays)  % history of variables
%      obs_history   (struct of arrays)  % history of observables
%      history_iterations (map)          % keep track of recorded
%                                          iterations
%
%
% Examples: 100 iterations of an arithmetic sequence of progression 1, 
%           with init value of 1.
%           Observables comprise the average of all terms.
%
%    model.variables.x = [1]; % all variables must be arrays
%                             % -> must encapsulate scalar values
%    model_solved = itmethod(model, arith_progess, 100, ...
%                            1, online_x_avg, 1, 1); 
%
%    function [model] = arith_progress(model)
%          model.variables.x = model.variables.x + 1;
%
%    function [avg] = online_x_avg(world)
%          if ~isfield(world.tmp.cumul_x)
%              world.tmp.cumul_x = zeros(size(world.variables.x);
%          end
%          world.tmp.cumul_x = world.observables.cumul + ...
%                              world.variables.x;
%          world.observables.avg_x = world.tmp.cumul / world.iterations;
% TODO: cut history to actual number of iterations
% TODO: manage recording of history iterations
% TODO: clarify iteration number for initialization
% TODO: test with different scenarios -> espacially with different
%       optional args
% TODO: add a pace to control computation of observables
    %% Argument processing
    ip_ = inputParser;
    ip_.addRequired('model', @isstruct);
    ip_.addRequired('update_variables', @isfunction);
    ip_.addRequired('max_iterations', @isnumeric);
    ip_.addOptional('var_history_pace', 0, @isnumeric);
    ip_.addOptional('compute_observables', @isfunction);
    ip_.addOptional('obs_cpt_start', 1, @isnumeric);
    ip_.addOptional('obs_history_pace', 0, @isnumeric);
    ip_.addOptional('stop_criterion', @(w) 0, @isfunction);
    ip_.addOptional('not_to_track', {}, @iscell);
    ip_.parse(model, update_variables, max_iterations, varargin{:});
    model = ip_.Results.model;
    update_vars = ip_.Results.update_variables;
    max_iterations = ip_.Results.max_iterations;
    stop_criterion = ip_.Results.stop_criterion;
    not_to_track = ip_.Results.not_to_track;
    vars_to_track = ~ismember(fieldnames(model.variables),  not_to_track);
    
    compute_observables = ip_.Results.compute_observables;
    obs_cpt_start = ip_.Results.obs_cpt_start;
    model.var_history_pace = ip_.Results.var_history_pace;
    model.obs_history_pace = ip_.Results.obs_history_pace;
    
    if ~exist('compute_observables','var')
       compute_observables_flag = false;% use a flag rather than testing if
                                        % compute_observables is defined
                                        % later because exist function 
                                        % is costly
    else
        compute_observables_flag = true;
    end
    
    %% Initializations and allocations
    if model.var_history_pace > 0
        max_size = max_iterations / model.var_history_pace;
        model.vars_history = allocate_history(model.variables, max_size, vars_to_track);
    end
    model.cur_var_hist_idx = 1;
    model.cur_obs_hist_idx = 1;
    model.iteration = 1;
    model.obs_iteration = 1; 
    
    model.var_hist_iterations = [];
    model = record_current_vars(model, vars_to_track);
    
    %% Main loop
    while model.iteration < max_iterations && ~stop_criterion(model)
        % Compute observables before updating to possibly take into account
        % initialization of variables
        if compute_observables_flag && model.iteration >= obs_cpt_start
            model = compute_observables(model);
            model.obs_iteration = model.obs_iteration + 1;
            if ~isfield(model, 'obs_history')
                % allocate history for observables 
                max_size = max_iterations / model.obs_history_pace - obs_cpt_start;
                obs_to_track = ~ismember(fieldnames(model.observables),  not_to_track);
                model.obs_history = allocate_history(model.observables,...
                                                     max_size, obs_to_track);
                model.obs_hist_iterations = [];
            end
            model = record_current_obs(model, obs_to_track);
        end
        model = update_vars(model);
        model.iteration = model.iteration + 1;
        model = record_current_vars(model, vars_to_track);
    end
    model.max_iterations = max_iterations;
    iterated_model = model;
end

function [history] = allocate_history(s, max_size, to_track)
    % Allocate history arrays for all fields in given struct *s*
    %
    % Return a struct with the same fields as *s* where each entry
    % is an array of size (*max_size*, shape of array in s)
    
    field_names = fieldnames(s);
    field_names = field_names(to_track); 
    history = struct();
    for i=1:numel(field_names)
        a = s.(field_names{i});
        history.(field_names{i}) = zeros([max_size, size(a)]);
    end
end

function [world] = record_current_vars(world, to_track)
    world = record(world, 'variables', 'vars_history', ...
                   'var_history_pace', 'cur_var_hist_idx', ...
                   'var_hist_iterations', to_track);
end

function [world] = record_current_obs(world, to_track)
    world = record(world, 'observables', 'obs_history', ...
                   'obs_history_pace', 'cur_obs_hist_idx', ...
                   'obs_hist_iterations', to_track);
end

function [model] = record(model, var_name, var_hist_name, pace_name, hist_idx_name, hist_its_name, to_track)
    % TODO: use Map to associate var_name with hist_cur_idx
    % -> lighten input args
    pace = model.(pace_name);
    if mod(model.iteration, pace) == 0 || ...
            model.iteration == 1 % always record start point  
        field_names = fieldnames(model.(var_name));
        field_names = field_names(to_track);
        for i=1:numel(field_names)
            if 0
                vh = model.(var_hist_name).(field_names{i});
                v = model.(var_name).(field_names{i});
                S.type = '()';
                S.subs{1} = model.(hist_idx_name);
                for j=1:ndims(v)
                    S.subs{j+1} = ':';
                end
                model.(var_hist_name).(field_names{i}) = subsasgn(vh, S, v);
            else
                switch ndims((model.(var_hist_name).(field_names{i})))
                    case 2
                        model.(var_hist_name).(field_names{i})(model.(hist_idx_name),:) = model.(var_name).(field_names{i});
                    case 3
                        model.(var_hist_name).(field_names{i})(model.(hist_idx_name),:,:) = model.(var_name).(field_names{i});
                    otherwise
                        error(['unhandled number of dims for history tracking of ' field_names{i}]);
                end
            end
%             disp([field_names{i} ':']);
%             disp(size(model.(var_hist_name).(field_names{i})));
        end
        model.(hist_its_name) = [model.(hist_its_name) model.iteration];
        model.(hist_idx_name) = model.(hist_idx_name) + 1;
    end
end