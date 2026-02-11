function [f,d,h] = ukernel_gauss(x,z,h,w)

    n_x = length(x);
    n_z = length(z);
    f   = zeros(length(x),1);
    d   = zeros(length(x),1);

%   Calculate the bandwidth
    if h <= 0
        s   = sqrt(var(z));
        q1  = prctile(z,25);
        q3  = prctile(z,75);
        iqr = q3 - q1;
        h   = 0.9*min([s, (iqr/1.34)])/n_z^0.2;
    end
    
    for i = 1:1:n_x
        arg = (x(i,1) - z)/h;
        kff = normpdf(arg,0,1);
        kfd = arg.*normpdf(arg,0,1); 
        f(i,1) = mean(kff.*w)/h;
        d(i,1) = mean(kfd.*w)/(h^2);
    end

end