function vprime = get_derivatives(model_data, ii)
% function vprime = get_derivatives runs the model code to compute
% derivates (e.g. time derivative of h_s) and secondary variables (e.g.
% phi_s).
% Inputs:
%   * outputs: model output structure returned by run_model.m
%   * ii: index of output to take derivatives of
% Returns:
%   vprime: structure with fields
%     vprime.phis_x:        x-component of gradient of sheet potential
%     vprime.phis_y:        y-component of gradient of sheet potential
%     vprime.phis_bndry:    sheet potential on element-edge boundaries
%     vprime.phic:          Channel potential
%
%     vprime.dhsd:          Time derivative of sheet thickness hs
%     vprime.dhcdt:         Time derivative of channel thickness hc
%     vprime.dHcdt:         Time derivative of channel incision Hc
%     vprime.dhsdt_div:     Divergence term component of dhsdt
%     vprime.dhcdt_div:     Divergence term component of dhcdt
%
%     vprime.qx_edge:       x-component of sheet flux on element-edge
%                            boundaries
%     vprime.qy_edge:       y-component of sheet flux on element-edge
%                            boundaries
%     vprime.qx_sheet:      x-component of sheet flux on element centroids
%     vprime.qy_sheet:      y-component of sheet flux on element centroids
%     vprime.qc:            Channel flux
%     vprime.exchange_frac: Mass-exchange fraction

% Construct state vector
outputs = model_data.outputs;
params = model_data.params;

v = [outputs.hs(:,ii); outputs.zs(:,ii);
     outputs.hc(:,ii); outputs.Hc(:,ii); outputs.phic(:,ii)];

% Now compute derivates
% if strcmp(params.channel_model, 'nonlinear')
vprime = rhs_sheet_and_channel(outputs.tt(ii), v, params.dmesh,...
            params, 'derivatives');
% elseif strcmp(params.channel_model, 'linear')
%     vprime = rhs_sheet_and_channel_linear(outputs.tt(ii), v, params.dmesh,...
%                 params, 'derivatives');
% end
