outs = load('outputs/baseline_more_moulins.mat');

mouts = load('outputs/spinup_refined_02_more_moulins.mat');
dmesh = outs.params.dmesh;

figure
edge_plot(gca, dmesh, abs(outs.outputs.qc(:, end)), palettes('-green-1'), 'none')

ii_moulin = load('data/moulins_M1.txt');
% ii_moulin = find(mouts.params.moulins);
% for ii=ii_moulin
plot(dmesh.tri.nodes(ii_moulin, 1), dmesh.tri.nodes(ii_moulin, 2), 'ro')
% end
% 
% z_nodes = shmip_elevation(dmesh.tri.nodes);
% 
% z_nodes(dmesh.tri.bmark>0) = Inf;   % No moulins on boundaries
% [Y, E] = discretize(z_nodes, 0:250:1500);
% 
% N_moulins = 0;
% moulin_density = [0.6, 0.45, 0.45, 0.2, 0.2, 0.2];
% % moulin_counts = [75, 90, 75, 50, 25, 15];
% ii_moulin = [];

% for binindex=1:length(E)-1
%     zindices = find(Y==binindex);
%     num_moulin = round(moulin_density(binindex)*length(zindices));
%     N_moulins = N_moulins + num_moulin;
%     disp('Number of nodes in elevation band:')
%     disp(length(zindices))
%     disp('Number of moulins in elevation band:')
%     disp(num_moulin)
%     
%     rng(2*binindex+4);
%     i_moulin = randi(length(zindices), 1, num_moulin);
%     ii_moulin = [ii_moulin, zindices(i_moulin)'];
% %     figure
% %     histogram(z_nodes(zindices))
% end
% N_moulins
% plot(dmesh.tri.nodes(ii_moulin, 1), dmesh.tri.nodes(ii_moulin, 2), 'bo')
% 