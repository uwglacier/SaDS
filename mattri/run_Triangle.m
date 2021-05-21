function run_Triangle(filename, triarea, degree_constraint, verbosity)
% run_Triangle(filename, triarea, degree_constraint, verbosity)
%
% Low level function.
%
% Runs Triangle which will write files in the directory of filename.
%
% If there is an .area file present, triarea will not be used.
%
% Input:
% - filename: Two possibilities:
%             - full filename of a .node or .poly file
%             - filename root.  If that is specified then a .area file needs to be present.
%               Used for refinements.

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if ~exist('triarea', 'var') || isempty(triarea)
    triarea = [];
end

if ~exist('degree_constraint', 'var') || isempty(degree_constraint)
    degree_constraint = 30;
end

if ~exist('verbosity', 'var') || isempty(verbosity)
    verbosity = 2;
end


current_dir = pwd;

%dir_of_triangle_executable = [getenv('HOME') , '/glads/numerics/matlab/mesh/triangle/'];
mfiledir = fileparts(mfilename('fullpath'));
if isunix()
    triangle_executable = [mfiledir, '/triangle/triangle'];
else
    triangle_executable = [mfiledir, '\triangle\triangle.exe'];
end

[path,filenameroot,extension] = fileparts(filename);
stdopts = '-CVDevj';
switch extension
  case '.node'
    command = [triangle_executable, ' ', stdopts ,' -q', num2str(degree_constraint), ' -a', num2str(triarea), ' ', filename];
  case '.poly'
    command = [triangle_executable, ' ', stdopts ,' -p -q', num2str(degree_constraint), ' -a', num2str(triarea), ' ', filename];
  case '' % with a .area file
    command = [triangle_executable, ' ', stdopts ,' -p -r -q', num2str(degree_constraint), ' -a ', filename ];
  otherwise
    if isnumeric(str2num(extension)) % we're refining meshes with .1, .2, etc in the filename root
        command = [triangle_executable, ' ', stdopts ,' -p -r -q', num2str(degree_constraint), ' -a ', filename ];
    else
        error(['Wrong input file extenstion: ',filename])
    end
end

if verbosity>0
    disp('Running Triangle with command:')
    disp(command)
end

[status,result] = unix(command);
if verbosity>1
    disp(result)
end
if strfind(result, 'No such file or directory')
    error('Cannot find triangle executable in glads_matlab/mesh/triangle/.  Have you built it?  See top of glads_matlab/mesh/triangle/README.')
end
if status~=0
    disp(result)
    error('Error running triangle')
end


% used flags
% -q minimum angle
% -a maximal area of triangle
% -D make conforming Delauny: - this ensures that the Voronoi is
% good.
% -V show detailed stats
% -e produce edge file
% -v produce voronoi files
% -j Jettisons vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).
%
% -p Triangulates a Planar Straight Line Graph (.poly file).
%
% -r refine

