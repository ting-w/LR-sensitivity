function [tSim, f, z,f_int,z_int, fz_int, count] = GT(model, tFinal, Ntr)
nu = model.nu;  %state change vector
network_dim = model.network_size;  %problem size
m = network_dim(1);
n = network_dim(2);
l = network_dim(3);
x = repmat(model.x0, 1, Ntr);  %initial population.
%z_old = zeros(Ntr, m);
f = model.obser(x);
z = zeros(l, Ntr);  %weight process
%obser = repmat(x, m, 1);
f_int = zeros(1, Ntr);
z_int = zeros(l, Ntr);
fz_int = zeros(l, Ntr); %m sensitivities
c = model.param;  %rate constant
tSim = zeros(1, Ntr);  %accumulated model time
finished = false(1, Ntr);
count = 0;
while any(~finished)
    count = count + 1;
    a = model.prop(x, c);
    jacob = model.jacob(x, c);
    ratio = model.ratio(x, c);
    a0 = cumsum(a);
    %    if a0(end, :) ~= 0  %there are reactions to fire
    rnd = rand(2, Ntr);
    tau = (1./a0(end, :)).*log(1./rnd(1, :));
    threshold = a0(end, :).*rnd(2, :);
    reacts = 1 + sum(a0 < repmat(threshold, m, 1));
    delx = nu(:, reacts);
    %update tau, size 1x100000
    tau(~finished) = min(tau(~finished),tFinal - tSim(~finished));
    %update time
    tSim(~finished) = tSim(~finished) + tau(~finished);
    z_old = z;
    %update the continuous part of weight process
    %reshape to make 3D to 2D: 1x100000x2 --> 100000x2
    z(:,~finished) = z(:,~finished) - reshape(sum(jacob(:,~finished,:),1),sum(~finished),l)'.*repmat(tau(~finished),l,1);
    f_int(:, ~finished) = f_int(:, ~finished) + f(:, ~finished).*tau(~finished);
    z_int(:, ~finished) = z_int(:, ~finished) + 0.5*(z_old(:, ~finished) + z(:, ~finished)).*repmat(tau(~finished),l,1);
    fz_int(:, ~finished) = fz_int(:, ~finished) + 0.5*repmat(f(:, ~finished),l,1).*repmat(tau(~finished),l,1).*(z_old(:,~finished)+z(:,~finished));  
    %update finished samples/trajectories
    %DO NOT update the discrete part of z before checking 'tSim >= tFinal'
    finished = (tSim >= tFinal);
    if any(~finished)
        effective_reacts = reacts(~finished);
        effective_ratio = ratio(:,~finished,:);
        %linear indexing: for every m elements, pick the right
        %effective_ratio in order to update Z_1, ..., Z_l
        index = (0 : l*(length(effective_reacts)) - 1)*m + repmat(effective_reacts, 1, l);
        %update the discrete part of weight process
        z(:,~finished) = z(:,~finished) + reshape(effective_ratio(index),length(effective_reacts),l)';
        %update state process
        x(:,~finished) = x(:,~finished) + delx(:, ~finished);
        f(:,~finished) = model.obser(x(:,~finished)); 
    else
        break;
    end
end