function out_struct = save_model_outputs(filename, params, Y)
% function save_model_outputs saves model outputs and parameters to a
% matlab structure for easy analysis later. Args are
%   * filename: where to save outputs to
%   * params: parameter structure used to run the model
%   * Y: Outputs straight from the ODE solver plus fields returned by
%        'derivatives' mode of solver
% 
% The resulting file has fields
% 
%   params: Structure listing all model parameters
%   outputs: Parameter providing model outputs. Struct with fields:
%                hs: [n_elements × n_outputs double]
%                zs: [n_elements × n_outputs double]
%                hc: [n_edges × n_outputs double]
%                Hc: [n_edges × n_outputs double]
%                zc: [n_edges × n_outputs double]
%                tt: [1 × n_outputs double]
%                Vm: [2 × n_outputs double]
%            phis_x: [n_elements × n_outputs double]
%            phis_y: [n_elements × n_outputs double]
%        phis_bndry: [n_edges × n_outputs double]
%             dhsdt: [n_elements × n_outputs double]
%             dhcdt: [n_edges × n_outputs double]
%             dHcdt: [n_edges × n_outputs double]
%         dhsdt_div: [n_elements × n_outputs double]
%         dhcdt_div: [n_edges × n_outputs double]
%           qx_edge: [n_edges × n_outputs double]
%           qy_edge: [n_edges × n_outputs double]
%          qx_sheet: [n_elements × n_outputs double]
%          qy_sheet: [n_elements × n_outputs double]
%                qc: [n_edges × n_outputs double]
%              phic: [n_edges × n_outputs double]
%          dphic_ds: [n_edges × n_outputs double]
%              Xi_c: [n_edges × n_outputs double]
%     exchange_frac: [n_edges × n_outputs double]
%          m_moulin: [n_moulins × n_outputs double]

% Initialize output structure
out_struct = struct;
% Put parameters into output structure
out_struct.params = params;

% Unpack model outputs
nouts=size(Y,1);
model_output.hs=zeros(params.dmesh.tri.n_elements,nouts);
model_output.zs=zeros(params.dmesh.tri.n_elements,nouts);
model_output.hc=zeros(params.dmesh.tri.n_edges,nouts);
model_output.Hc=zeros(params.dmesh.tri.n_edges,nouts);
model_output.phic=zeros(params.dmesh.tri.n_edges,nouts);
model_output.tt = params.tt;

for ii=1:nouts
    [hs,zs,hc,Hc,zc,vm]=unpack_state_vector(params.dmesh,Y(ii,:)');
    model_output.hs(:,ii)=hs;
    model_output.zs(:,ii)=zs;
    model_output.hc(:,ii)=hc;
    model_output.Hc(:,ii)=Hc;
    model_output.phic(:,ii)=zc;
    model_output.Vm(:,ii)=vm;
   
    % Compute derivatives
    deriv_load.outputs = model_output;
    deriv_load.params = params;
    dx = get_derivatives(deriv_load, ii);
    for jj=1:length(params.output_fields)
        field = params.output_fields{jj};
        if isfield(model_output, field)
            model_output.(field) = [model_output.(field), dx.(field)];
        else
            model_output.(field) = [dx.(field)];
        end            
    end
end

out_struct.outputs=model_output;
out_struct.outputs.tt = out_struct.params.tt;

if ~strcmp('none', filename)
     try
        [fpath, f2, f3] = fileparts(filename);
        % Check if the folder exists
        if ~isfolder(fpath)
            mkdir(fpath)
        end

        if ~isfile(filename) || params.overwrite
%            if isfile(filename)
%                delete(filename)
%                disp(sprintf('Deleted file %s', filename))
%            end
%            curdir = pwd;
%            cd(fpath);
            save(filename,'-struct','out_struct');
%            disp('Attemped to save file')
%            cd(curdir);
        else
            d = datetime;
            d.Format = '_uuuu-MM-dd_mm_ss';
            split_name = split(filename, '.');
            split_name{1} = [split_name{1}, char(d)];
            new_filename = join(split_name, '.');
            save(new_filename{1}, '-struct', 'out_struct');
        end
     catch
         d = datetime;
         d.Format = '_uuuu-MM-dd_mm_ss';
         new_filename = ['outputs', char(d), '.mat'];
         save(new_filename, '-struct', 'out_struct');
         warning('Could not save outputs as %s, saved as %s', filename, new_filename)
     end
end
