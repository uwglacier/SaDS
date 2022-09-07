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

write_nc(outs, ncfile)
