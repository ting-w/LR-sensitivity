clear;
model = sixD_system;
l = model.network_size(3);  % parameter #
Ntr = 1000;    % sample # in each batch
Batch = 100;   % total batch #
t_final = 1.0e+4;
for itr = 1 : 100
    [t_arr{itr}, f_arr{itr}, z_arr{itr}, f_int_arr{itr}, z_int_arr{itr},...
        fz_int_arr{itr}, count{itr}] = GT(model, t_final, Ntr);
end

%save 'GT_sixD_1e4';
% f = mean(f_arr, 2);
% z = mean(z_arr, 2);
% %ergodic average for x
% f_int = mean(f_int_arr/t_final, 2);
% %ergodic average for z
% z_int = mean(z_int_arr/t_final, 2);
% %GT array
% GT_arr = repmat(f_int_arr/t_final, l, 1).*z_arr;
% %CGT array
% CGT_arr = repmat(f_int_arr/t_final-f_int,l,1).*z_arr;
% %WGT array
% WGT_arr = fz_int_arr/t_final;
% %WCGT array
% %WCGT_arr = WGT_arr - repmat(x_int_arr/t_final, 2, 1).*(z_int_arr'/t_final);
% WCGT_arr = WGT_arr - repmat(f_int,l,Ntr).*(z_int_arr/t_final);
% 
% %GT estimate
% GT = mean(GT_arr,2);
% %CGT estimate
% CGT = mean(CGT_arr,2);
% %CGT = GT - x_int*z';
% %weighted estimate
% WGT = mean(WGT_arr, 2);
% %weighted centered estimator
% WCGT = mean(WCGT_arr, 2);
% %WCGT = WGT-x_int*z_int';
% 
% %var of GT
% Var_GT = sum((GT_arr - repmat(GT,1,Ntr)).^2,2)/(Ntr-1);
% %var of CGT
% Var_CGT = sum((CGT_arr - repmat(CGT,1,Ntr)).^2,2)/(Ntr-1);
% %var of WGT
% Var_WGT = sum((WGT_arr - repmat(WGT,1,Ntr)).^2,2)/(Ntr-1);
% %var of WCGT
% Var_WCGT = sum((WCGT_arr - repmat(WCGT,1,Ntr)).^2,2)/(Ntr-1);
% 
