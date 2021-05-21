function mesh_plot_tri(ax, meshtri, bmarky, bmarky_nodes, label, color)
% mesh_plot_tri(ax, meshtri, bmarky, bmarky_nodes, label, color)
%
% Plots a triangular mesh in blue with the boundary edges in different
% color if option bmarky==1: black for bmark_edge==1 and green for
% bmark_edge==2.
%
% If label == 'nodes' then label nodes
%    label == 'edges' then label edges

% This m-file is part of the mattri meshing package.
% Copyright (c) 2014, Mauro A Werder
% All rights reserved.
% Licensed under a BSD 2-Clause License, see LICENCE file

if isempty(ax)
    figure
    ax = gca;
end
    
if ~exist('bmarky', 'var') || isempty(bmarky)
    bmarky = 0;
end
if ~exist('bmarky_nodes', 'var') || isempty(bmarky_nodes)
    bmarky_nodes = 0;
end
if ~exist('label', 'var') || isempty(label)
    label = 'none';
end
if ~exist('color', 'var') || isempty(color)
    color = 'b';
end

if bmarky<2
axes(ax);
triplot(meshtri.connect, meshtri.nodes(:,1), meshtri.nodes(:,2), color);
hold on;
plot(meshtri.nodes(:,1), meshtri.nodes(:,2), [color,'.']);
end

cols = ['g', 'r', 'c', 'm', 'y', 'k']; %, 'w'];
% make boundary different colors
if bmarky>0
    % edges with bmark odd
    bm1 = logical(mod(meshtri.bmark_edge,2));
    eg1 = get_edges(bm1);
    plot(eg1(:,1),eg1(:,2), 'g')%, plot_args{:})

% $$$     bm2 = meshtri.bmark_edge==2;
% $$$     eg1 = get_edges(bm2);   
% $$$     plot(eg1(:,1),eg1(:,2), 'b' , plot_args{:})
 
    % edges with bmark even and > 0
    for ii = 1:length(cols)*10
        colind = mod(ii,length(cols))+1;
        bm2 = meshtri.bmark_edge==(ii)*2;

        eg1 = get_edges(bm2);   
        plot(eg1(:,1),eg1(:,2), cols(colind))% , plot_args{:})
    end

% $$$     bm2 = ~logical(mod(meshtri.bmark_edge,2)) & meshtri.bmark_edge>2;
% $$$     %    bm2 = meshtri.bmark_edge==4;
% $$$     eg1 = get_edges(bm2);   
% $$$     plot(eg1(:,1),eg1(:,2), 'k', plot_args{:})
end

if bmarky_nodes>0
    % nodes
    bm2 = meshtri.bmark==2;

    plot(meshtri.nodes(bm2,1),meshtri.nodes(bm2,2), 'g.')

    bm1 = logical(mod(meshtri.bmark,2));  % odd ones
    plot(meshtri.nodes(bm1,1),meshtri.nodes(bm1,2), 'k.', 'markersize',20)
    plot(meshtri.nodes(bm1,1),meshtri.nodes(bm1,2), 'r.', 'markersize',10)    


    % nodes with bmark even and > 0
    for ii = 1:length(cols)*2
        colind = mod(ii,length(cols))+1;
        bm2 = meshtri.bmark==(ii)*2;
        plot(meshtri.nodes(bm2,1),meshtri.nodes(bm2,2), [cols(colind),'.'], 'markersize',15)   
    end
end


% labels
switch lower(label)
  case 'nodes'
    labs = cellfun(@num2str,num2cell(1:size(meshtri.nodes,1)),'UniformOutput', false);
    text(meshtri.nodes(:,1), meshtri.nodes(:,2), labs);
  case 'edges'
    labs = cellfun(@num2str,num2cell(1:size(meshtri.connect_edge,1)),'UniformOutput', false);
    text(meshtri.edge_midpoints(:,1), meshtri.edge_midpoints(:,2), labs);
end


hold off;
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edges = get_edges(bmarks)

    eg1a = meshtri.nodes(meshtri.connect_edge(bmarks,1),:);
    eg1b = meshtri.nodes(meshtri.connect_edge(bmarks,2),:);   
    eg1 = zeros(size(eg1a,1)*3,2)*NaN; % we need NaN inserted
                                       % between edges to stop
                                       % wrong edges being drawn
    eg1(1:3:end,:) = eg1a;
    eg1(2:3:end,:) = eg1b ;   
    edges = eg1;
    end
end