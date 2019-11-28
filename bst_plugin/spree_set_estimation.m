function model = nirs_model_set_estimation(model, vars_to_estimate)
% Enable the estimation of the given variable and disable the estimation
% of all others
%
% Inputs ([]s are optional)
%     (struct) model
%        main data container with the following fields:
%           - (struct of arrays) variables
%             Model variables to be updated. Given *variable_name*
%             must be a field of this struct.
%     (str) vars_to_name
%        Variables to estimate. Must be fields in model.variables.  


var_names = fieldnames(model.variables);


for iv=1:length(var_names)
    model.variables_to_estimate.(var_names{iv}) = 0;
end

if isfield(model, 'cst_to_estimate')
    cst_names = fieldnames(model.cst_to_estimate);
    for ic=1:length(cst_names)
        model.cst_to_estimate.(cst_names{ic}) = 0;
    end
end

for iv=1:length(vars_to_estimate)
    vname = vars_to_estimate{iv};
    if ismember(vname, var_names)
        model.variables_to_estimate.(vname) = 1;
    elseif isfield(model, 'cst_to_estimate') && ismember(vname, cst_names)
        model.cst_to_estimate.(vname) = 1;
    else
        throw(MException('NIRSDynaError:WrongVariableName', ...
                         ['Cannot find variable ' vname]));
    end
end

end