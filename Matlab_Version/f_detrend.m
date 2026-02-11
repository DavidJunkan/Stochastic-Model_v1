function [yld_dt,yld_t,yld_dv] = f_detrend(yld,method)

    n_crop = width(yld);
    yld_t  = zeros(size(yld));
    yld_c  = zeros(size(yld));
    yld_dv = zeros(size(yld));
    yld_dt = zeros(size(yld));

    for i = 1:1:n_crop
        if method == 1
            v_max = max(yld(:,i));
            v_min = min(yld(:,i));
            yld_t(:,i) = linspace(v_min,v_max,length(yld))';
            yld_c(:,i) = yld(:,i) - yld_t(:,i);
        elseif method == 2
            T = (1:1:length(yld))';
            X = [ones(length(yld),1),T];
            reg_ts = fitlm(X, yld(:,i), 'Intercept', false);
            yld_t(:,i) = X*reg_ts.Coefficients.Estimate;
            yld_c(:,i) = yld(:,i) - yld_t(:,i);
        elseif method == 3
            [y_t,y_c]  = hpfilter(yld(:,i));
            yld_t(:,i) = y_t;
            yld_c(:,i) = y_c;
        end
        yld_dv(:,i) = yld_c(:,i)./yld_t(:,i);
        yld_dt(:,i) = (yld(:,i)./yld_t(:,i))*yld_t(length(yld),i);
    end

end
