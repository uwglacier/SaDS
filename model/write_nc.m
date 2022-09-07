function write_nc(S, ncfile)
% mat2nc: Savemodel outputs in nc format
%
% mat2nc(S, ncfile) saves model outputs in struct S to file ncfile and
% metadata XML file with the same name

% Add extension
[f1, f2, ~] = fileparts(ncfile);
metafile = [f1, f2, '.meta'];

%% Write nc file
% Find dimensions
n_elements = length(S.outputs.hs(:, 1));
n_edges = length(S.outputs.hc(:, 1));
n_moulins = length(S.outputs.m_moulin(:, 1));
n_times = length(S.outputs.tt);

% Start creating netcdf variables
nccreate(ncfile, 'hs', 'Dimensions', {'n_elements', n_elements, 'n_times', n_times})
ncwrite(ncfile, 'hs', S.outputs.hs)

nccreate(ncfile, 'zs', 'Dimensions', {'n_elements', n_elements, 'n_times', n_times})
ncwrite(ncfile, 'zs', S.outputs.zs)

nccreate(ncfile, 'hc', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
ncwrite(ncfile, 'hc', S.outputs.hs)

nccreate(ncfile, 'Hc', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
ncwrite(ncfile, 'Hc', S.outputs.Hc)

nccreate(ncfile, 'phic', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
ncwrite(ncfile, 'phic', S.outputs.phic)

nccreate(ncfile, 'tt', 'Dimensions', {'n_times', n_times})
ncwrite(ncfile, 'tt', S.outputs.tt)

nccreate(ncfile, 'm_moulin', 'Dimensions', {'n_moulin', n_moulins, 'n_times', n_times})
ncwrite(ncfile, 'm_moulin', S.outputs.m_moulin)

nccreate(ncfile, 'qc', 'Dimensions', {'n_edges', n_edges, 'n_times', n_times})
ncwrite(ncfile, 'qc', S.outputs.qc)

%% Write META file

% Need to change datatypes of some of the fields
print_struct = S.params;
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
