function p_bins = gen_pbin(p_simu,bins,idx_bin)

%   If idx_m == 1, use one bin for all crops, else if idx_m == 0, using
%   different bins for different crops
    if idx_bin == 1
        bin_s   = bins(:,1);
        width_s = width(p_simu);
        bin_mat = kron(ones(1,width_s),bin_s);
    else
        width_s = min([width(p_simu),width(bins)]);
        bin_mat = bins(:,1:width_s);
        p_simu  = p_simu(:,1:width_s);
    end

    p_bins = zeros(size(bin_mat));
    prob_s = (1:1:length(p_simu))'/length(p_simu);

    for i = 1:1:width_s
        [in_val,~]  = sort(p_simu(:,i), 'ascend');
        cum_p_bins  = cumsum(bin_mat(:,i));
        p_bins(:,i) = interp1(prob_s,in_val,cum_p_bins);
    end

end