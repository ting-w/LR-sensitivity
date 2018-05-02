function mymodel = model
mymodel.nu = [1    -1     0     0     0     0     0     0     0     0     0     0;
              0     0     1    -1    -2     2     0     0     0     0    -1     1;
              0     0     0     0     1    -1     0     0     0     0     0     0;
              0     0     0     0     0     0     1    -1     0     0     0     0;
              0     0     0     0     0     0     0     0     1    -1    -1     1;
              0     0     0     0     0     0     0     0     0     0     1    -1];
%initial population
mymodel.x0 = [0; 0; 0; 0; 0; 0];
%rate constant
%[k_r, phi, k_dr, k_p, k_dp, k_1, k_2, k_3, k_4]
mymodel.param = [1 60 0.1 1 0.5 0.02 0.08 0.02 0.1];
%propensity function
mymodel.prop = @(X, param) propFun(X, param);
%Jacobian
mymodel.jacob = @(X, param) jacobFun(X, param);
%ratio
mymodel.ratio = @(X, param) ratio(X, param);
mymodel.obser = @(X) obser(X);
mymodel.network_size = [12; 6; 9]; %[m,n,l] reaction types #; species #; paramter #
end

%vectorized propensity functions
function output = propFun(X, param)
output = [param(1)*param(2)^4 ./ (param(2)^4 + X(3,:).^4);
    param(3)*X(1, :);
    param(4)*X(1, :);
    param(5)*X(2, :);
    param(6)*X(2, :).*(X(2, :) - 1);
    param(7)*X(3, :);
    param(1)*param(2)^2 ./ (param(2)^2 + X(6, :).^2);
    param(3)*X(4, :);
    param(4)*X(1, :);
    param(5)*X(5, :);
    param(8)*X(2, :).*X(5, :);
    param(9)*X(6, :)];
end


%We only consider sensitivity wrt x1
function output = jacobFun(X, param)
h = length(X(1,:));
output(:, :, 1) = [param(2)^4 ./ (param(2)^4 + X(3,:).^4); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h); 
    param(2)^2 ./ (param(2)^2 + X(6, :).^2); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 2) = [4*param(1)*param(2)^3./(param(2)^4 + X(3, :).^4) - 4*param(1)*param(2)^7./(param(2)^4 + X(3, :).^4).^2; zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    2*param(1)*param(2)./(param(2)^2 + X(6, :).^2) - 2*param(1)*param(2)^3./(param(2)^2 + X(6, :).^2).^2; zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 3) = [zeros(1, h); X(1, :); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); X(4, :); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 4) = [zeros(1, h); zeros(1, h); X(1, :); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); X(1, :);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 5) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    X(2, :); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    X(5, :); zeros(1, h); zeros(1, h)];

output(:, :, 6) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); X(2, :).*(X(2, :) - 1); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 7) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); X(3, :);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 8) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); X(2, :).*X(5, :); zeros(1, h)];

output(:, :, 9) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); X(6, :)];
end


function output = ratio(X, param)
h = length(X(1,:));

output(:, :, 1) = [ones(1, h)/param(1); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h); 
    ones(1, h)/param(1); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 2) = [4/param(2) - 4*param(2)^3./(param(2)^4 + X(3, :).^4); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    2/param(2) - 2*param(2)./(param(2)^2 + X(6, :).^2); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 3) = [zeros(1, h); ones(1, h)/param(3); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); ones(1, h)/param(3); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 4) = [zeros(1, h); zeros(1, h); ones(1, h)/param(4); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); ones(1, h)/param(4);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 5) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    ones(1, h)/param(5); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    ones(1, h)/param(5); zeros(1, h); zeros(1, h)];

output(:, :, 6) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); ones(1, h)/param(6); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 7) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); ones(1, h)/param(7);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h)];

output(:, :, 8) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); ones(1, h)/param(8); zeros(1, h)];

output(:, :, 9) = [zeros(1, h); zeros(1, h); zeros(1, h); 
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); zeros(1, h);
    zeros(1, h); zeros(1, h); ones(1, h)/param(9)];

%row:reaction types
%column: sample numbers
%page:rate constant
end

function output = obser(X)
output = X(6,:);
end
