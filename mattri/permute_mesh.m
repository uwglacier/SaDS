function mesh = permute_mesh(mesh, per, per_edge)
%  mesh = permute_mesh(mesh, per, per_edge)
%
% Low level function
% 
% Applies a permutation of the node numbering and a different one to the edge numbering.
%
% mesh.tri.connect -- list nodes for elements
% mesh.tri.connect_edge -- list nodes for edges
% mesh.tri.connect_edge_inv -- list edges for nodes
% mesh.tri.connect_edge_el -- lists elements for edges
%
% mesh.vor.connect -- lists edges for cells
% mesh.vor.connect_edge -- lists cells for edges
% mesh.vor.connect_cell -- lists cells for cells

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file


%%%%%%%%
%% nodes

% inverse permutation http://blogs.mathworks.com/loren/2007/08/21/reversal-of-a-sort/
tmp = 1:mesh.tri.n_nodes;
inv_per(per) = tmp;

%% tri part
% nodes
mesh.tri.nodes = mesh.tri.nodes(per,:);
mesh.tri.area_nodes = mesh.tri.area_nodes(per);
mesh.tri.bmark = mesh.tri.bmark(per);

% elements
mesh.tri.connect = inv_per(mesh.tri.connect); %[inv_per(mesh.tri.connect(:,1))', inv_per(mesh.tri.connect(:,2))', inv_per(mesh.tri.connect(:,3))'];
% permutation of elements: none as just the "name" of the nodes changes
mesh.tri.area = mesh.tri.area; % same
if isfield(mesh, 'vor')
    mesh.tri.connect_edge_el = mesh.tri.connect_edge_el; % same as edge and element numbers stay
end

% edges
mesh.tri.connect_edge = inv_per(mesh.tri.connect_edge); %[inv_per(mesh.tri.connect_edge(:,1))', inv_per(mesh.tri.connect_edge(:,2))'];

% permutation of edges: none as just the "name" of the nodes changes
mesh.tri.bmark_edge = mesh.tri.bmark_edge;  % same
mesh.tri.edge_length = mesh.tri.edge_length; %same
mesh.tri.edge_midpoints = mesh.tri.edge_midpoints; % same


%% vor part
if isfield(mesh, 'vor')
    % cells correpond to nodes and thus need to be reorderd
    mesh.vor.bmark_cell = mesh.vor.bmark_cell(per);
    mesh.vor.connect = mesh.vor.connect(per); % connection between cells & edges
    
    % this is the trickest:
    mesh.vor.connect_cell = mesh.vor.connect_cell(per);
    for ii = 1:length(mesh.vor.connect_cell)
        mesh.vor.connect_cell{ii} = inv_per(mesh.vor.connect_cell{ii});
    end
end
    

if nargin==4 && ~isempty(per_edge)
    %%%%%%%%
    %% edges
    
    tmp = 1:mesh.tri.n_edges;
    inv_per_edge(per_edge) = tmp;
    
    %% tri
    mesh.tri.bmark_edge = mesh.tri.bmark_edge(per_edge);
    mesh.tri.edge_length = mesh.tri.edge_length(per_edge);
    mesh.tri.edge_midpoints = mesh.tri.edge_midpoints(per_edge);
    
    mesh.tri.connect_edge = mesh.tri.connect_edge(per_edge,:);
    mesh.tri.connect_edge_el = mesh.tri.connect_edge_el(per_edge,:);
    
    if isfield(mesh, 'vor')
        %% vor
        % voronoi edges correspond to tri-edges and thus need reordering
        %
        % There is some iffyness because of the trimming of the Voronoi mesh, which thus has more
        % edges than the triangular mesh!  The additional edges are added to the end thus we just
        % ignore those.
        
        mesh.vor.connect_edge(tmp,:) = mesh.vor.connect_edge(per_edge,:);
        mesh.vor.bmark_edge(tmp) = mesh.vor.bmark_edge(per_edge);
        for ii = 1:length(mesh.vor.connect)
            ind = mesh.vor.connect{ii}<=mesh.tri.n_edges;
            mesh.vor.connect{ii}(ind) = inv_per_edge(mesh.vor.connect{ii}(ind));
        end
    end
end

%%%%%%%%
% finishing up

% recalculate connect_edge_inv:
for nd =1:mesh.tri.n_nodes
    tmp1 = (sort(find(mesh.tri.connect_edge(:,1)==nd)))';
    tmp2 = (sort(find(mesh.tri.connect_edge(:,2)==nd)))';
    connect_edge_inv{nd} = [tmp1, tmp2];
end
mesh.tri.connect_edge_inv = connect_edge_inv;
