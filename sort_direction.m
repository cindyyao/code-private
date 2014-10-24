function ds_out = sort_direction(ds_in)

% ds_in: structure contain all ds information
% xyao 
% 2014-08-12

for j = 1:length(ds_in.rho)
    [theta_seq, I_seq] = sort(ds_in.theta{j}(1, :));
    r = ds_in.rho{j};
    R = ds_in.RHO{j};
    rho_seq = r(:, I_seq);
    RHO_seq = R(:, I_seq);
    ds_in.rho{j} = rho_seq;
    ds_in.RHO{j} = RHO_seq;
    ds_in.theta{j} = repmat(theta_seq, size(ds_in.theta{j}, 1), 1);
end
ds_out = ds_in;
end