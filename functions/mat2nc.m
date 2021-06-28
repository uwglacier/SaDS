function mat2nc(fpath, varargin)
% mat2nc: Convert .mat files to netCDF format
%
% mat2nc(fpath) converts file 'fpath.mat' to 'fpath.nc', additionally
% creating metadata file fpath.meta.
%
% met2nc(fpath, ncpath) saves the netcdf file as ncpath.nc and the metadata
% file ncpath.meta

% Read in .mat file
outs = load(fpath);

% Parse inputs for NC path
if isempty(varargin)
    ncpath = fpath;
else
    ncpath = varargin{1};
end

% Add extension
ncfile = [ncpath, '.nc'];
metafile = [ncpath, '.meta'];

%% Write nc file
% Find dimensions
n_elements = length(outs.outputs.hs(:, 1));
n_edges = length(outs.outputs.hc(:, 1));
n_moulins = length(outs.outputs.m_moulin(:, 1));
n_times = length(outs.outputs.tt);

% % Start creating netcdf variables
% nccreate(ncfile, 'hs', 'Dimensions', {'n_elements', n_elements, 'n_times', n_times})
% ncwrite(ncfile, 'hs', outs.outputs.hs)
% 
% nccreate(ncfile, 'zs', 'Dimensions', {'n_elements', n_elements, 'n_times', n_times})
% ncwrite(ncfile, 'zs', outs.outputs.zs)
% 
% nccreate(ncfile, 'hc', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
% ncwrite(ncfile, 'hc', outs.outputs.hs)
% 
% nccreate(ncfile, 'Hc', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
% ncwrite(ncfile, 'Hc', outs.outputs.Hc)
% 
% nccreate(ncfile, 'phic', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
% ncwrite(ncfile, 'phic', outs.outputs.phic)
% 
% nccreate(ncfile, 'tt', 'Dimensions', {'n_times', n_times})
% ncwrite(ncfile, 'tt', outs.outputs.tt)
% 
% nccreate(ncfile, 'm_moulin', 'Dimensions', {'n_moulin', n_edges, 'n_times', n_times})
% ncwrite(ncfile, 'm_moulin', outs.outputs.m_moulin)
% 
% nccreate(ncfile, 'qc', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
% ncwrite(ncfile, 'qc', outs.outputs.qc)

%% Write META file

% Need to change datatypes of some of the fields
print_struct = outs.params;
print_struct.moulins = find(full(print_struct.moulins)==1);
print_struct.q_neumann = [];
print_struct.qN = [];
print_struct.dmesh = [];
print_struct.tt = [];
print_struct.wc = [];
print_struct.solver_opts.OutputFcn = char(print_struct.solver_opts.OutputFcn);
print_struct.solver_opts.stats = [];

print_struct.ms = char(print_struct.ms);
print_struct.msc = char(print_struct.msc);
print_struct.mc = char(print_struct.mc);
print_struct.fc = char(print_struct.fc);

writestruct(print_struct, metafile, 'FileType', 'xml');
