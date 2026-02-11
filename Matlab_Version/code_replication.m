clear
clc

%% Basic Settings
n_simu   = 15000;
y_base   = 2024;
mu       = 0;
val_miss = -9.999;
n_spline = 1000;    

%% Read in the price data
%  In the price matrix, each column indicates one crop type, users can define their own orders, but please make sure the order of crop prices matches the order of crop yields
%  Current Order: Corn, Sorghum, Barley, Oats, Wheat, Rice, Up Cotton, Soybeans, Rice(LG), Rice(MG)

p_input  = 'price_input.xlsx';
T        = readmatrix(p_input, 'Sheet', 1, 'Range', 'A2:K13');
T(isnan(T)) = val_miss;
y_p_det  = T(:,1);
p_det    = T(:,2:width(T));
n_crop   = width(p_det);
n_yp_det = length(y_p_det);

T        = readmatrix(p_input, 'Sheet', 2, 'Range', 'B2:K3');
sigma    = T(1,:);
sigma(isnan(sigma)) = 0;
time_adj = T(2,:);
time_adj(isnan(time_adj)) = 1;

bins     = readmatrix('bins.xlsx', 'Sheet', 1, 'Range', 'B2:B35');
seeds    = readmatrix('seedmat3.csv');

T        = readmatrix(p_input, 'Sheet', 3, 'Range', 'A2:K13');
T(isnan(T)) = val_miss;
y_a_det  = T(:,1);
acre_det = T(:,2:width(T));

T        = readmatrix(p_input, 'Sheet', 4, 'Range', 'A2:K13');
T(isnan(T)) = val_miss;
y_q_det  = T(:,1);
prod_det = T(:,2:width(T));

T        = readmatrix(p_input, 'Sheet', 5, 'Range', 'A2:E26');
T(isnan(T)) = val_miss;
year_p   = T(:,1);
p_dv_dat = T(:,2:width(T));

%% Price Simulation
p_val    = p_det(1,:);
seed_val = seeds(1,1);
[p_simu_dv,p_simu] = gen_logn(mu,sigma,time_adj,p_val,n_simu,seed_val);
check_p  = mean(p_simu) - p_val; 

%% Read in the acre and yield historical data
%  Crop order in acre and production data should also be consistent with the order in the price input

n_year   = 50;
year_q   = (y_base - n_year + 1:y_base)';
acre_dat = zeros(n_year,n_crop);
prod_dat = zeros(n_year,n_crop);

for i = 1:1:n_crop
    opts   = detectImportOptions('acre_input.xlsx','Sheet',i);
    opts.VariableNamingRule = 'preserve';
    T_acre = readtable('acre_input.xlsx', opts);
    A_acre = T_acre{:,2:4};
    opts   = detectImportOptions('production_input.xlsx','Sheet',i);
    opts.VariableNamingRule = 'preserve';
    T_prod = readtable('production_input.xlsx', opts);
    A_prod = T_prod{:,2:4};
    for j = 1:1:n_year
        set_1 = find(A_acre(:,1) == year_q(j,1) & A_acre(:,2) == 99);
        set_2 = find(A_prod(:,1) == year_q(j,1) & A_prod(:,2) == 99);
        if ~isempty(set_1)
            acre_dat(j,i) = A_acre(set_1(1),3);
        else
            acre_dat(j,i) = val_miss;
        end
        if ~isempty(set_2)
            prod_dat(j,i) = A_prod(set_2(1),3);
        else
            prod_dat(j,i) = val_miss;
        end
    end
end

%% Yield Simulation, Kernel
yld_val     = prod_dat./acre_dat;
% dt_method = 2: Simple linear trend; dt_method = 2: OLS with linear trend;
% dt_method = 3: HP filter
dt_method   = 2;
[yld_dt, yld_t, yld_dv] = f_detrend(yld_val,dt_method);
yld_simu    = zeros(n_simu,n_crop);
yld_simu_dv = zeros(n_simu,n_crop);

for i = 1:1:n_crop
    yld_min   = min(yld_dt(:,i));
    yld_max   = max(yld_dt(:,i));
    ext_gap   = 2*std(yld_dt(:,i));
  % ext_gap   = 0;
    yld_ssp   = linspace(yld_min - ext_gap, yld_max + ext_gap, n_spline)';
  
   %[pdf,yld_ssp]   = kde(yld_dt(:,i), NumPoints = n_spline);
    [pdf,~,~] = ukernel_gauss(yld_ssp,yld_dt(:,i),0,1);
    P         = pdf.*(max(yld_ssp) - min(yld_ssp))/(n_spline - 1);
    P         = cumsum(P);
    P_l       = min(P);
    P_h       = max(P);
    U         = rand(n_simu,1);
    for j = 1:1:n_simu
        if U(j,1) <= P_l
            id_l = 1;
            id_h = 2;
        elseif U(j,1) >= P_h
            id_l = n_spline - 1;
            id_h = n_spline;
        else
            idx_P  = P - U(j,1);
            select = find(idx_P > 0);
            id_l   = select(1,1) - 1;
            id_h   = select(1,1);
        end
        yld_simu(j,i) = yld_ssp(id_l,1)*(P(id_h,1) - U(j,1))/(P(id_h,1) - P(id_l,1)) +  yld_ssp(id_h,1)*(U(j,1) - P(id_l,1))/(P(id_h,1) - P(id_l,1));
        yld_simu(j,i) = max(yld_simu(j,i),min(yld_ssp));
        yld_simu(j,i) = min(yld_simu(j,i),max(yld_ssp));
    end
    yld_simu_dv(:,i)  = (yld_simu(:,i) - mean(yld_simu(:,i)))./mean(yld_simu(:,i));
end

mean_real  = mean(yld_dt);
mean_simu  = mean(yld_simu);
std_real   = std(yld_dt);
std_simu   = std(yld_simu);

kurt_real  = kurtosis(yld_dt);
kurt_simu  = kurtosis(yld_simu);
skew_real  = skewness(yld_dt);
skew_simu  = skewness(yld_simu);

check_simu = [mean_real; std_real; kurt_real; skew_real; mean_simu; std_simu; kurt_simu; skew_simu];

%% Capula Method
year_cor   = intersect(year_p, year_q);
p_dv_cor   = zeros(length(year_cor),n_crop);
yld_dv_cor = zeros(length(year_cor),n_crop);

for i = 1:1:length(year_cor)
    id_p = find(year_p == year_cor(i,1));
    id_q = find(year_q == year_cor(i,1));
    p_dv_cor(i,:)   = p_dv_dat(id_p,:);
    yld_dv_cor(i,:) = yld_dv(id_q,:);
end

%  Note that the price deviation of rice, LG rice, and MG rice are identical
data_cor = [p_dv_cor(:,1:n_crop-2),yld_dv_cor];
rho      = corr(data_cor);
pp       = chol(rho);
XX_cor   = randn(n_simu,length(rho));
RXX_cor  = XX_cor*pp;
UXX_cor  = normcdf(RXX_cor,0,1);

simu_raw = [p_simu_dv(:,1:n_crop-2),yld_simu_dv];
simu_cor = zeros(size(simu_raw));
prob_ref = (1:1:n_simu)'/n_simu;

for i = 1:1:width(data_cor)
    x_val = sort(simu_raw(:,i),'ascend');
    for j = 1:1:n_simu
        if UXX_cor(j,i) <= 1/n_simu
            id_l = 1;
            id_h = 2;
        else
            idx_P  = prob_ref - UXX_cor(j,i);
            select = find(idx_P > 0);
            id_l   = select(1,1) - 1;
            id_h   = select(1,1);
        end
        simu_cor(j,i) = x_val(id_l,1)*(prob_ref(id_h,1) - UXX_cor(j,i))/(prob_ref(id_h,1) - prob_ref(id_l,1)) +  x_val(id_h,1)*(UXX_cor(j,i) - prob_ref(id_l,1))/(prob_ref(id_h,1) - prob_ref(id_l,1));
    end
    simu_cor(:,i) = simu_cor(:,i) - mean(simu_cor(:,i));
end

rho_simu = corr(simu_cor);
check_rho = rho - rho_simu;
p_simu_dv_cor   = simu_cor(:,1:n_crop - 2);
yld_simu_dv_cor = simu_cor(:,n_crop - 1:width(data_cor));

%% Now use the baseline determinstic values
%  You are free to decide which year as your target

y_target   = 2024;
p_target   = p_det(y_p_det == y_target,:);
arc_target = acre_det(y_a_det == y_target,:);
yld_target = prod_det(y_q_det == y_target,:)./arc_target;

yld_det_simu_iid  = zeros(n_simu,n_crop);
yld_det_simu_cor  = zeros(n_simu,n_crop);
prod_det_simu_iid = zeros(n_simu,n_crop);
prod_det_simu_cor = zeros(n_simu,n_crop);
p_det_simu_iid    = zeros(n_simu,n_crop);
p_det_simu_cor    = zeros(n_simu,n_crop);

for i = 1:1:n_crop
    yld_det_simu_iid(:,i)  = yld_simu_dv(:,i)*yld_target(1,i) + yld_target(1,i);
    yld_det_simu_cor(:,i)  = yld_simu_dv_cor(:,i)*yld_target(1,i) + yld_target(1,i);
    prod_det_simu_iid(:,i) = yld_det_simu_iid(:,i)*arc_target(1,i);
    prod_det_simu_cor(:,i) = yld_det_simu_cor(:,i)*arc_target(1,i);
    if i <= n_crop - 2
        id_p_dv = i;
    else
        id_p_dv = 6;    % i_crop = 6 is rice
    end
    p_det_simu_iid(:,i) = p_simu_dv(:,i)*p_target(1,i) + p_target(1,i);
    p_det_simu_cor(:,i) = p_simu_dv_cor(:,id_p_dv)*p_target(1,i) + p_target(1,i);
end

mean_simu_iid  = [mean(p_det_simu_iid), mean(yld_det_simu_iid)];
std_simu_iid   = [std(p_det_simu_iid),  std(yld_det_simu_iid)];
kurt_simu_iid  = [kurtosis(p_det_simu_iid), kurtosis(yld_det_simu_iid)];
skew_simu_iid  = [skewness(p_det_simu_iid), skewness(yld_det_simu_iid)];

mean_simu_cor  = [mean(p_det_simu_cor), mean(yld_det_simu_cor)];
std_simu_cor   = [std(p_det_simu_cor),  std(yld_det_simu_cor)];
kurt_simu_cor  = [kurtosis(p_det_simu_cor), kurtosis(yld_det_simu_cor)];
skew_simu_cor  = [skewness(p_det_simu_cor), skewness(yld_det_simu_cor)];

check_cor = [mean_simu_iid; std_simu_iid; kurt_simu_iid; skew_simu_iid;
             mean_simu_cor; std_simu_cor; kurt_simu_cor; skew_simu_cor];


%% Generate Distribution by Bins
%  same target year
idx_bin  = 1;
bins_p_cor   = gen_pbin(p_det_simu_cor,bins,idx_bin);
bins_yld_cor = gen_pbin(yld_det_simu_cor,bins,idx_bin);


%% Plot (1578)
figure
subplot(3,3,1)
crop_x1 = 1;
crop_x2 = 1 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
legend([rho_1 rho_2], {' Data',' Simulated'},'Location','best','FontSize',12)
xlabel('Corn Price','Fontsize',12)
ylabel('Corn Yield','Fontsize',12)


subplot(3,3,2)
crop_x1 = 5;
crop_x2 = 5 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Wheat Price','Fontsize',12)
ylabel('Wheat Yield','Fontsize',12)


subplot(3,3,3)
crop_x1 = 7;
crop_x2 = 7 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Cotton Price','Fontsize',12)
ylabel('Cotton Yield','Fontsize',12)


subplot(3,3,4)
crop_x1 = 8;
crop_x2 = 8 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('OrRd'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Soybean Price','Fontsize',12)
ylabel('Soybean Yield','Fontsize',12)



%  Yield Correlation
figure
subplot(3,3,1)
crop_x1 = 5 + 8;
crop_x2 = 1 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
%legend([rho_1 rho_2], {' Data',' Simulated'},'Location','best','FontSize',8)
xlabel('Corn Yield','Fontsize',12)
ylabel('Wheat Yield','Fontsize',12)


subplot(3,3,2)
crop_x1 = 7 + 8;
crop_x2 = 1 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Corn Yield','Fontsize',12)
ylabel('Cotton Yield','Fontsize',12)


subplot(3,3,3)
crop_x1 = 8 + 8;
crop_x2 = 1 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Corn Yield','Fontsize',12)
ylabel('Soybean Yield','Fontsize',12)


subplot(3,3,4)
crop_x1 = 7 + 8;
crop_x2 = 5 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Wheat Yield','Fontsize',12)
ylabel('Cotton Yield','Fontsize',12)


subplot(3,3,5)
crop_x1 = 8 + 8;
crop_x2 = 5 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Wheat Yield','Fontsize',12)
ylabel('Soybean Yield','Fontsize',12)

subplot(3,3,6)
crop_x1 = 8 + 8;
crop_x2 = 7 + 8;
xx      = [simu_cor(:,crop_x2), simu_cor(:,crop_x1)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
xlabel('Cotton Yield','Fontsize',12)
ylabel('Soybean Yield','Fontsize',12)



%   Price Correlation
figure
subplot(3,3,1)
crop_x1 = 5;
crop_x2 = 1;
xx      = [simu_cor(:,crop_x1), simu_cor(:,crop_x2)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
%legend([rho_1 rho_2], {' Data',' Simulated'},'Location','best','FontSize',8)
ylabel('Corn Price','Fontsize',12)
xlabel('Wheat Price','Fontsize',12)


subplot(3,3,2)
crop_x1 = 7;
crop_x2 = 1;
xx      = [simu_cor(:,crop_x1), simu_cor(:,crop_x2)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
ylabel('Corn Price','Fontsize',12)
xlabel('Cotton Price','Fontsize',12)


subplot(3,3,3)
crop_x1 = 8;
crop_x2 = 1;
xx      = [simu_cor(:,crop_x1), simu_cor(:,crop_x2)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
ylabel('Corn Price','Fontsize',12)
xlabel('Soybean Price','Fontsize',12)


subplot(3,3,4)
crop_x1 = 7;
crop_x2 = 5;
xx      = [simu_cor(:,crop_x1), simu_cor(:,crop_x2)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
ylabel('Wheat Price','Fontsize',12)
xlabel('Cotton Price','Fontsize',12)


subplot(3,3,5)
crop_x1 = 8;
crop_x2 = 5;
xx      = [simu_cor(:,crop_x1), simu_cor(:,crop_x2)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('YlGnBu'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
ylabel('Wheat Price','Fontsize',12)
xlabel('Soybean Price','Fontsize',12)

subplot(3,3,6)
crop_x1 = 8;
crop_x2 = 7;
xx      = [simu_cor(:,crop_x1), simu_cor(:,crop_x2)];
n_bin   = 20;
min_x   = min(min(xx(:,1)));
max_x   = max(max(xx(:,1)));
min_y   = min(min(xx(:,2)));
max_y   = max(max(xx(:,2)));
width_x = (max_x - min_x)/(n_bin-1);
width_y = (max_y - min_y)/(n_bin-1);
Xedges  = [1.1*min_x min_x:width_x:max_x 1.1*max_x];
Yedges  = [1.1*min_y min_y:width_y:max_y 1.1*max_y];
h = histogram2(xx(:,1),xx(:,2),Xedges,Yedges,'EdgeAlpha',0.2,...
    'Normalization','probability',...
    'DisplayStyle','tile','ShowEmptyBins','on');
colormap(slanCM('PuBuGn'))
caxis([0 0.03]);
colorbar('eastoutside')
zz      = linspace(0.9*min_x,0.9*max_x,100);
yy_data = zz*rho(crop_x1,crop_x2);
yy_simu = zz*rho_simu(crop_x1,crop_x2);
hold on
rho_1 = plot(zz,yy_data,'k-','Linewidth', 2);
rho_2 = plot(zz,yy_simu,'r--','Linewidth', 2);
hold off
ylabel('Cotton Price','Fontsize',12)
xlabel('Soybean Price','Fontsize',12)