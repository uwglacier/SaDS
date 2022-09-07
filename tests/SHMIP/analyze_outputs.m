% Script to analyze SaDS model outputs
% This script uses the functions in the /functions/ directory

set_paths

% First we read in the outputs
outs = load('./outputs/baseline.mat');
dmesh = outs.params.dmesh;

% Plot sheet water depth
figure('Units', 'centimeters', 'Position', [4, 4, 20, 8])
element_plot(dmesh, outs.outputs.hs(:, end), 'EdgeColor', 'none')
axis image

% Now make it a nicer plot
xlabel('Easting (m)')
ylabel('Northing (m)')
cmocean('dense')
cb = colorbar;
cb.Label.String = 'h_s (m)';

print('fig_hs', '-dpng', '-r600')

% Plot channel flow
figure('Units', 'centimeters', 'Position', [3, 3, 20, 8])
edge_plot(gca, dmesh, abs(outs.outputs.qc(:, end)), palettes('-green-1'), [0, 15])

% Now make it a nicer plot
xlabel('Easting (m)')
ylabel('Northing (m)')
cb = get(gca, 'ColorBar');
cb.Label.String = 'h_s (m)';

print('fig_qc', '-dpng', '-r600')
