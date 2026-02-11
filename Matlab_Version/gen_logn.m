function [p_dev,p_simu] = gen_logn(mu,sigma,time_adj,p_val,n_simu,seed)

%   Check if the dimension match
    width_p  = width(p_val);
    width_s  = width(sigma);
    width_t  = width(time_adj);
    width_v  = min([width_p, width_s, width_t]);
    p_val    = p_val(:,1:width_v);
    time_adj = time_adj(:,1:width_v);
    sigma    = sigma(:,1:width_v);

%   Simulate the prices
    sigmaT   = sqrt(time_adj).*sigma;
    mean_ZT  = (mu - (sigma.^2)/2).*time_adj;
    ln_p     = log(p_val);
    ln_p_mat = kron(ones(n_simu,1),ln_p);
    mean_mat = kron(ones(n_simu,1),mean_ZT);
    rng(seed);
    rv_norm  = randn(n_simu, width(p_val));
    ln_pr    = ln_p_mat + mean_mat + rv_norm.*sigmaT;
    p_simu   = exp(ln_pr);
    p_dev    = (p_simu - exp(ln_p_mat))./exp(ln_p_mat);
    diff_m   = p_val - mean(p_simu);
    p_simu   = p_simu + kron(ones(n_simu,1),diff_m);

end