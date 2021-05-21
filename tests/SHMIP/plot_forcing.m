% Script to plot the SHMIP forcing. Plots a Hovmoller-style plot of the
% forcing for the entire elevation range, the total melt at each elevation,
% and the maximum melt to show the diurnal cycle.

tt = (110:0.1:256)*86400;
z = 0:10:1500;

[X, Y] = meshgrid(tt, z);

melt = shmip_melt(Y, X, false);

figure('Units', 'inches', 'Position', [2, 2, 8, 5])

axes('Position', [0.1, 0.25, 0.7, 0.7])
set(gca, 'FontSize', 9)
hold on
imagesc(tt/86400, z, 86400*melt)
% contour(X/86400, Y, 86400*melt, [1e-3, 0.01, 0.1], 'k')
% shading flat
cmocean('rain')
cbar = colorbar('Location', 'northoutside');
caxis([0, 0.12])
xlim([110, 256])
ylim([0, 1500])
set(gca, 'YTick', [0, 500, 1000, 1500])
set(gca, 'XTickLabel', {'', '', '', '', '', '', '', ''})
set(gca, 'XTick', [110, 130, 150, 170, 190, 210, 230, 250])
% colorbar_label(cbar, 'Melt')
text(0.1 + (0.5*0.75), 1.25, 'Mean melt (m w.e./day)', 'Units', 'normalized', 'HorizontalAlignment', 'center')
% xlabel('Day of year')
ylabel('Elevation (m asl.)', 'FontSize', 9)
text(0.025, 0.975, 'a', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9)

pos = get(gca, 'Position');
cm = get(gca, 'Colormap');

axes('Position', [0.8, 0.25, 0.15, pos(4)])
set(gca, 'FontSize', 9)
plot(sum(melt*(tt(2) - tt(1)), 2), z, 'Color', cm(170, :))
set(gca, 'YTick', [0, 500, 1000, 1500])
set(gca, 'YTickLabel', {})
xlabel('Total melt (m w.e.)', 'FontSize', 9)
set(gca, 'XTick', [2, 4, 6, 8, 10])
text(0.05, 0.975, 'b', 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9)
set(gca, 'Box', 'off')

axes('Position', [0.1, 0.1, 0.7, 0.15])
set(gca, 'FontSize', 9)
melt_diurnal = shmip_melt(0, tt, true);
plot(tt/86400, 86400*melt_diurnal, 'Color', cm(170, :))
set(gca, 'YDir', 'reverse')
ylim([0, 0.24])
xlim([110, 256])
set(gca, 'Box', 'off')
set(gca, 'YTick', [0.1, 0.2])
text(0.025, 0.9, 'c', 'Units', 'normalized', 'VerticalAlignment', 'top')
xlabel('Day of year', 'FontSize', 9)
ylabel('Melt (m w.e./day)', 'FontSize', 9)
set(gca, 'XTick', [110, 130, 150, 170, 190, 210, 230, 250])

print('melt_parameterization', '-dpng', '-r600')
print('melt_parameterization', '-dpdf', '-r600', '-opengl')
