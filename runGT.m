clear;
model = model;
l = model.network_size(3);  % parameter #
Ntr = 1000;    % sample # in each batch
Batch = 100;   % total batch #
t_final = 1.0e+4;
for itr = 1 : 100
    [t_arr{itr}, f_arr{itr}, z_arr{itr}, f_int_arr{itr}, z_int_arr{itr},...
        fz_int_arr{itr}, count{itr}] = GT(model, t_final, Ntr);
end

