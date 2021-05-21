function combined_plot(ax,plot_data)
% combined_plot. Overlay edge variable onto element variable.
%
% combined_plot(ax, plot_data) makes the plot on axes ax. plot_data is a
% structure with fields
%             dmesh: [1×1 struct]
%     Cdata_element: [n_elements×1 double]
%      cmap_element: [N×3 double]
%        Cdata_edge: [n_edges×1 double]
%        caxis_edge: 'none' or [cmin, cmax]
%          varargin: Cell array passed to edge_plot
%         cmap_edge: [N×3 double]
%          edge_min: double
% See also edge_plot, element_plot

dmesh = plot_data.dmesh;

% Plot elements by colouring their faces
patch('Faces', dmesh.tri.connect, 'Vertices', dmesh.tri.nodes,...
        'FaceVertexCData', plot_data.Cdata_element, 'FaceColor', 'flat', 'EdgeColor', 'none');
colormap(ax, plot_data.cmap_element);
cb = colorbar;
% cb.Label.String = plot_data.cb_label_element;
% Put colorbar on north outside of the plot
cb.Location = 'northoutside';
axis image

% Make another axes for the edges
ax2 = axes;
% Hide the second axes
ax2.Visible = 'off';
% Colour edges with heavier lines than usual
plot_data.Cdata_edge(plot_data.Cdata_edge<plot_data.edge_min) = nan;
edge_plot(ax2,plot_data.dmesh,plot_data.Cdata_edge,plot_data.cmap_edge,plot_data.caxis_edge,'Lineweights', [1, 3], 'vmin', plot_data.edge_min)
axis image
% cb2 = get(ancestor(ax2, 'axes'), 'ColorBar');
% cb2.Label.String = plot_data.cb_label_edge;
% cb2.Label.Rotation = 0;
% cb2.Label.HorizontalAlignment = 'left';

% Make axes have same limits
linkaxes([ax, ax2])

% Make axes have identical positions
axpos = get(ax, 'Position');
ax2pos = get(ax2, 'Position');
set([ax, ax2], 'Position', [axpos(1), axpos(2), ax2pos(3), axpos(4)]);
