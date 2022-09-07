function colorbar_label(cb, labelstr)
% colorbar_label. Label a colorbar.
%
% colorbar_label(cb, labelstr) labels the Colorbar instance cb with the
% string labelstr. This function only works for colorbars oriented
% vertically in their default position. Alternatively, label the colorbar
% directly:
%   cb = colorbar;
%   cb.Label.String = 'foo';

cb.Label.String = labelstr;
cb.Label.Units = 'normalized';
cb.Label.Position = [0.5, 1.25, 0];
cb.Label.HorizontalAlignment = 'center';
cb.Label.Rotation = 0;
end
