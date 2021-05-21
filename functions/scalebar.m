function scalebar(ax, majorticks, majorticklabels, minorticks, minorticklabels)
% scalebar. Draw a GIS-style scalebar on the specified axes
%
% scalebar(ax, majorticks, majorticklabels, minorticks, minorticklabels)
% Add scalebar to axes ax. Place major ticks at positions majorticks, and
% label them with labels majorticklabels. Same for minorticks and
% minorticklabels. Can set to [] and {} to not plot major or minor
% ticks.
%
% EXAMPLE: Scalebar with no minor ticks
%
% figure
% element_plot(dmesh, hs(:, end))
% set(gca, 'Visible', 'off') % No axes labels so may want a scale bar
% scalebar(gca, 1e3*[0, 5, 10], {'0', '5 km', '10 km'})
%
% See also element_plot, edge_plot

majorticklength = 0.025;
minorticklength = 0.015;
lw = 0.67;
minorlw = 0.5;

xlims = get(gca, 'Xlim');
ylims = get(gca, 'Ylim');
ymin = ylims(1);
dy = ylims(2) - ylims(1);
xmin = xlims(1);

hold(ax, 'on')
plot(ax, [xmin, xmin + majorticks(end)], [ymin, ymin], 'k', 'LineWidth', 0.67)

for ii=1:length(majorticks)
    abslength = majorticklength*dy;
    plot([xmin + majorticks(ii), xmin + majorticks(ii)], [ymin, ymin + abslength], 'k', 'LineWidth', lw)
    if ~isempty(majorticklabels)
        text(xmin + majorticks(ii), ymin + 1.25*abslength, majorticklabels{ii},...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline')
    end
end

for ii=1:length(minorticks)
    abslength = minorticklength*dy;
    plot([xmin + minorticks(ii), xmin + minorticks(ii)], [ymin, ymin + abslength], 'k', 'LineWidth', minorlw)
    if ~isempty(minorticklabels)
        text(xmin + minorticks(ii), ymin + abslength, minorticklabels{ii},...
            'HorizontalAlignment', 'center')
    end
end
