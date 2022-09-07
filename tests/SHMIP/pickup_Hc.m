function Hc = pickup_Hc(pk_file)
outs = load(pk_file);
Hc = outs.outputs.Hc(:, end);
Hc(Hc<0.25) = 0.25; % Simulate snow melt
