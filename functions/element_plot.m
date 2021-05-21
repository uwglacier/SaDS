function element_plot(dmesh,Cdata, varargin)
% element_plot. Plot values stored on elements.
%
% element_plot(dmesh, Cdata, varargin). Mesh object dmesh, Cdata is the
% data to plot (N_elements x 1), and varargin arguments are passed to
% Matlab's patch function
%
% element_plot is really just a thin wrapper to the patch function with
% an easier call signature
%
% See also patch, edge_plot

patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
        'FaceVertexCData', Cdata, 'FaceColor', 'flat', varargin{:});
