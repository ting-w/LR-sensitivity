clear;
% load('GT_sixD_5e3.mat');
% GT{1} = GT;
% WGT{1} = WGT;
% CGT{1} = CGT;
% WCGT{1} = WCGT;
% Var_GT{1} = Var_GT;
% Var_WGT{1} = Var_WGT;
% Var_CGT{1} = Var_CGT;
% Var_WCGT{1} = Var_WCGT;
% 
% load('GT_sixD_1e4.mat');
% GT{2} = GT;
% WGT{2} = WGT;
% CGT{2} = CGT;
% WCGT{2} = WCGT;
% Var_GT{2} = Var_GT;
% Var_WGT{2} = Var_WGT;
% Var_CGT{2} = Var_CGT;
% Var_WCGT{2} = Var_WCGT;
% 
% load('GT_sixD_1_5e4.mat');
% GT{3} = GT;
% WGT{3} = WGT;
% CGT{3} = CGT;
% WCGT{3} = WCGT;
% Var_GT{3} = Var_GT;
% Var_WGT{3} = Var_WGT;
% Var_CGT{3} = Var_CGT;
% Var_WCGT{3} = Var_WCGT;
% 
% load('GT_sixD_2e4.mat');
% GT{4} = GT;
% WGT{4} = WGT;
% CGT{4} = CGT;
% WCGT{4} = WCGT;
% Var_GT{4} = Var_GT;
% Var_WGT{4} = Var_WGT;
% Var_CGT{4} = Var_CGT;
% Var_WCGT{4} = Var_WCGT;
% 
% load('GT_sixD_2_5e4.mat');
% GT{5} = GT;
% WGT{5} = WGT;
% CGT{5} = CGT;
% WCGT{5} = WCGT;
% Var_GT{5} = Var_GT;
% Var_WGT{5} = Var_WGT;
% Var_CGT{5} = Var_CGT;
% Var_WCGT{5} = Var_WCGT;


load('GT_sixD_2_5e4.mat');
%sample size
Ntr = Batch * Ntr;

%concatenate each batch to one single array
count = cell2mat(count);
f_arr = cell2mat(f_arr);
f_int_arr = cell2mat(f_int_arr);
fz_int_arr = cell2mat(fz_int_arr);
z_arr = cell2mat(z_arr);
z_int_arr = cell2mat(z_int_arr);


f = mean(f_arr, 2);
z = mean(z_arr, 2);
%ergodic average for x
f_int = mean(f_int_arr/t_final, 2);
%ergodic average for z
z_int = mean(z_int_arr/t_final, 2);
%GT array
GT_arr = repmat(f_int_arr/t_final, l, 1).*z_arr;
%CGT array
CGT_arr = repmat(f_int_arr/t_final-f_int,l,1).*z_arr;
%WGT array
WGT_arr = fz_int_arr/t_final;
%WCGT array
%WCGT_arr = WGT_arr - repmat(x_int_arr/t_final, 2, 1).*(z_int_arr'/t_final);
WCGT_arr = WGT_arr - repmat(f_int,l,Ntr).*(z_int_arr/t_final);

%GT estimate
GT = mean(GT_arr,2);
%CGT estimate
CGT = mean(CGT_arr,2);
%CGT = GT - x_int*z';
%weighted estimate
WGT = mean(WGT_arr, 2);
%weighted centered estimator
WCGT = mean(WCGT_arr, 2);
%WCGT = WGT-x_int*z_int';

%var of GT
Var_GT = sum((GT_arr - repmat(GT,1,Ntr)).^2,2)/(Ntr-1);
%var of CGT
Var_CGT = sum((CGT_arr - repmat(CGT,1,Ntr)).^2,2)/(Ntr-1);
%var of WGT
Var_WGT = sum((WGT_arr - repmat(WGT,1,Ntr)).^2,2)/(Ntr-1);
%var of WCGT
Var_WCGT = sum((WCGT_arr - repmat(WCGT,1,Ntr)).^2,2)/(Ntr-1);

