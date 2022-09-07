% Script computes a hydrograph for a synthetic DEM using flow routing with
% a constant velocity

%% Inputs
% Flow velocity
v = 0.02;

dmesh = load('../meshes/phys_mesh.mat');

% Synthetic elevation map
mx = 0.005;
my = 0.0025;

% Node elevation
z = mx*dmesh.tri.nodes(:,1) + my*abs(dmesh.tri.nodes(:,2));

% Element elevation
z_el = mx*dmesh.tri.elements(:,1) + my*abs(dmesh.tri.elements(:,2));

% Melt rates (constant 5 cm/day at term., 2.5 cm/day at top) for half a day
m=@(zz,t) 0.05/86400 * (1 - 0.5*(zz - min(z))./(max(z) - min(z))) * heav(t).*heav(86400/2 - t);

%% Compute travel time T
[L, paths] = mesh_flow_routing(dmesh, z);
T = L./v;

%% Plot travel time
figure('Position',[680 558 600 300])
hold on
patch('Faces',dmesh.tri.connect,'Vertices',dmesh.tri.nodes,'FaceVertexCData',T/3600,'FaceColor','flat');
xlabel('x (m)')
ylabel('y (m)')
axis image
title('Flow routing travel time')
colorbar
colormap(gca,cmocean('tempo'))
text(1.07,1.12,'T (hr)','units','normalized')
print('flow_routing_travel_time','-dpng','-r400')

%% Compute hydrograph
dt = 600;
t = 0;
tend = 2*86400;
tt = t:dt:tend;
R = zeros(tend/dt + 1,1);
jj = 1;
while t <= tend
    for ii=1:dmesh.tri.n_elements
        melt = m(z_el(ii), t - T(ii));
        R(jj) = R(jj) + melt*dmesh.tri.area(ii);
    end
    jj = jj + 1;
    t = t + dt;
end

figure
plot(tt/3600, R)
title('Moulin discharge')
xlabel('Hours')
ylabel('Discharge (m^3/s)')
print('flow_routing_discharge','-dpng','-r400')