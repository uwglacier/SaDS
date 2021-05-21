function status = validate_params(params)
% validate_params validates the parameter structure before each simulation
%
% status = validate_params(params) checks if each field in params is valid
% (e.g. defined by the model), and that all required parameters are set.
%
% status = 1 means the check passed, and status = 0 means the check failed.
% Each invalid option is displayed as a warning.
%
% See also get_default_params

% Pass unless we find a problem
status = 1;

% Fields the user must specify that are NOT set by get_default_params
req_fields = {'tt', 'ms', 'msc', 'mc', 'fc'};

% Find the fields that are set by get_default_params
def_fields = fields(get_default_params(params.dmesh));

input_fields = fields(params);

%% Global checks ----------------------------------------------------------
% First, check that each field in params is recognized
for ii=1:length(input_fields)
    input_field = input_fields{ii};
    is_default = ismember(input_field, def_fields);
    is_reqd =  ismember(input_field, req_fields);
    if ~is_default && ~is_reqd
        warning('Input field `%s` not recognized', input_field);
        status = 0;
    end
end

% Second, check that each required field is given
for jj=1:length(req_fields)
    req_field = req_fields{jj};
    if ~ismember(req_field, input_fields)
        warning('Required field %s not provided', req_field)
        status = 0;
    end
end

%% Now check some cases individually --------------------------------------

% Ensure timestep is provided if using odeRK solver
if strcmp(params.solver, 'odeRK')
    if ~isfield(params.solver_opts, 'dt')
        warning('Timestep not specified with odeRK solver')
        status = 0;
    end
end

end
