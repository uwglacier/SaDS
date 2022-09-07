function xy = points_in_tri(nodes, varargin)
% points_in_tri returns a list of coordinates inside the given triangle
%
% USAGE:
% xy = points_in_tri(nodes)
%
% xy = points_in_tri(nodes, h) uses horizontal resolution h (default 100)
%

% Unpack node coordinates
x1 = nodes(1,1); x2 = nodes(2,1); x3 = nodes(3,1);
y1 = nodes(1,2); y2 = nodes(2,2); y3 = nodes(3,2);

% Generate mesh
xmin = min(nodes(:,1));
xmax = max(nodes(:,1));
ymin = min(nodes(:,2));
ymax = max(nodes(:,2));

if isempty(varargin)
    h = 100;
else
    h = varargin{1};
end

xs = xmin-h:h:xmax+h;
ys = ymin-h:h:ymax+h;

[xx, yy] = meshgrid(xs, ys);

in_tri = 0*xx;

for mm=1:size(xx, 1)
    for nn=1:size(xx,2)
        x = xx(mm, nn);
        y = yy(mm, nn);
        
        % Barycentric coordinate method
        a = ((y2-y3)*(x-x3) + (x3-x2)*(y-y3))/((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3));
        b = ((y3-y1)*(x-x3) + (x1-x3)*(y-y3))/((y2-y3)*(x1-x3) + (x3-x2)*(y1-y3));
        c = 1 - a - b;
        
        if 0<=a && 1>=a && 0<=b && 1>=b && 0<=c && 1>=c
            in_tri(mm, nn) = 1;
        end
    end
end

xtri = xx(in_tri==1);
ytri = yy(in_tri==1);
xy = [xtri, ytri];
