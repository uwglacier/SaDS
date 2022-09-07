% Compare the difference between gradient and absolute routing methods

moulins_abs = load('optimized_moulins_baseline_02.txt');
moulins_grad= load('gradient_optimized_moulins_baseline_02.txt');

dmesh = load('../SaDS/meshes/shmip_refined_mesh_02.mat');

outs = load('../SaDS/cases/shmip/sensitivity/outputs/baseline_opt.mat');

figure
edge_plot(gca, dmesh, abs(outs.outputs.qc(:, end)), palettes('-green-1'), 'none')
hold on
% for ii=1:dmesh.tri.n_elements
%     nodes = dmesh.tri.connect(ii, :);
%     nodes = [nodes, nodes(1)];
%     plot(dmesh.tri.nodes(nodes, 1), dmesh.tri.nodes(nodes, 2), 'Color', [0.5, 0.5, 0.5])
% end
% axis image

commons = intersect(moulins_abs, moulins_grad);
plot(dmesh.tri.nodes(moulins_abs, 1), dmesh.tri.nodes(moulins_abs, 2), 'ro')
plot(dmesh.tri.nodes(moulins_grad,1), dmesh.tri.nodes(moulins_grad,2), 'bo')
plot(dmesh.tri.nodes(commons, 1), dmesh.tri.nodes(commons, 2), 'go')
