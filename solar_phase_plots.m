close all; clear; clc;

% The least squares solution will take n number of measurements (at least
% 3) to compute the position of the spacecraft

% Define constant parameters
moon_angle = 20; % incidence angle with moon surface (degrees)
radius_Moon = 1.7374e6; % m
altitude0 = 35000;
altitudefinal = 97000;
delta_altitude = 10000; % How much to iterate altitude by

solar_vec = [20 50 80 110 140];

for i = 1:length(solar_vec)
moon_angle = solar_vec(i);
altitude_sc = altitude0;
fprintf('%d\n', i)
z = 1;
while altitude_sc <= altitudefinal


distance_sc = altitude_sc+radius_Moon; % m

% Calculate position vector of spacecraft used to locate craters using some
% random unit vector
r_sc_I = [0.4811; 0.2580; 0.8379]*(distance_sc); % m

% Calculate matrix of mean and std
[angular_errors,repeat_matrix_detections,x_2, y_2, z_2] = angular_error_calc(altitude_sc, moon_angle);

% Max index value
rep_matrix_size = size(repeat_matrix_detections);
num_craters = rep_matrix_size(1);
fprintf('Max index: %d\n', num_craters);
is_pass = false; % used to ensure index is chosen within the size of the array

% Display list of craters
% for j = 1:rep_matrix_size(1)
%     fprintf('%d -- sc_inc = %.2f -- crater_radius = %.2f\n', j, repeat_matrix_detections(j,3),repeat_matrix_detections(j,4))
% end

desired_indices = 1:num_craters; % list of indices to use for crater generation

% Extract crater quantities associated with desired indeces
std_devs = repeat_matrix_detections(desired_indices,2); % rad
means = repeat_matrix_detections(desired_indices,1); % rad
sc_bearings = repeat_matrix_detections(desired_indices, 3); % rad
crater_radius = repeat_matrix_detections(desired_indices, 4); % m

idx_Nan = find(isnan(std_devs)); % finds indices of Nan values

% Remove Nan values
std_devs(idx_Nan) = []; 
means(idx_Nan) = [];
sc_bearings(idx_Nan) = [];
crater_radius(idx_Nan) = [];

idx_zero = find(std_devs == 0); % finds indices of 0 std value

% Remove zero values
std_devs(idx_zero) = []; 
means(idx_zero) = [];
sc_bearings(idx_zero) = [];
crater_radius(idx_zero) = [];

% Number of times to run simulation

N = 10000; % to use for Monte Carlo

% Pre allocate space for errors
% mean_LS = zeros(N, 1);    - Used for first n craters
% std_LS = zeros(N, 1);
% 
% mean_WLS = zeros(N, 1);
% std_WLS = zeros(N, 1);

% for j = 1:jend

% Need to put above implementation in crater_pos, organize and just
% extracxt i number of craters

% The azimuth and elevation angles are not the same for each call of the
% function, thats why it needs to be added within the function, just add it
% there, sort, and extract the best i craters, use i as input to function

% Pre allocate space
% LS_errors = zeros(num_craters, N); % For n best craters
% WLS_errors = zeros(num_craters, N); 
LS_errors = zeros(1, N); % For altitude/convergence
WLS_errors = zeros(1, N);
mean_LS_con = zeros(1, N);
std_LS_con = zeros(1, N);
mean_WLS_con = zeros(1, N);
std_WLS_con = zeros(1, N);


% for i = 1:N % Monte Carlo
% for i = 1 % need minimum three measurements
    for j = 1:N % num iterations
    % Pass statistic values 
    [BCWLS_e, WLS_e, LS_e, r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearings, std_devs, means,crater_radius);
    
        LS_errors(:,j) = LS_e(end); % end only for altitude variation
        WLS_errors(:,j) = WLS_e(end);
        BCWLS_errors(:,j) = BCWLS_e(end);

        % crater_LS(j,:) = LS_e;
        % crater_WLS(j,:) = WLS_e;
        % crater_BCWLS(j,:) = BCWLS_e;

        % mean_LS_con(j) = mean(LS_errors(end,1:j),2);
        % std_LS_con(j) = std(LS_errors(end,1:j),0,2);
        % 
        % mean_WLS_con(j) = mean(WLS_errors(end,1:j),2);
        % std_WLS_con(j) = std(WLS_errors(end,1:j),0,2);
    
    end
mean_LS(i,z) = mean(LS_errors,2);
std_LS(i,z) = std(LS_errors,0,2);
% mean_crat_LS = mean(crater_LS,1);
% std_crat_LS = std(crater_LS,0,1);

mean_WLS(i,z) = mean(WLS_errors,2);
std_WLS(i,z) = std(WLS_errors,0,2);
% mean_crat_WLS = mean(crater_WLS,1);
% std_crat_WLS = std(crater_WLS,0,1);

mean_BCWLS(i,z) = mean(BCWLS_errors,2);
std_BCWLS(i,z) = std(BCWLS_errors,0,2);
% mean_crat_BCWLS = mean(crater_BCWLS,1);
% std_crat_BCWLS = std(crater_BCWLS,0,1);

fprintf("Iteration: %.d\n", altitude_sc)

% Update iterators
altitude_sc = altitude_sc + delta_altitude; % add 500 m each iteration
z = z + 1;
% end


% % Scale errors by altitude
% scaled_LS(:,j) = LS_errors(:,j)/altitude_sc * 100; 
% scaled_WLS(:,j) = WLS_errors(:,j)/altitude_sc * 100;
end

end



% Plot histograms for the scaled errors

% figure
% subplot(2,2,1)
% hold on
% histogram(scaled_LS(:,1), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5)
% histogram(scaled_WLS(:,1), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5, 'LineStyle','--')
% title('Scaled Least Squares Errors')
% xlabel('Error (%)')
% ylabel('Probability Density')
% 
% subplot(2,2,2)
% hold on
% histogram(scaled_LS(:,2), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5)
% histogram(scaled_WLS(:,2), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5, 'LineStyle','--')
% title('Scaled Least Squares Errors')
% xlabel('Error (%)')
% ylabel('Probability Density')
% 
% subplot(2,2,3)
% hold on
% histogram(scaled_LS(:,3), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5)
% histogram(scaled_WLS(:,3), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5, 'LineStyle','--')
% title('Scaled Least Squares Errors')
% xlabel('Error (%)')
% ylabel('Probability Density')
% 
% subplot(2,2,4)
% hold on
% histogram(scaled_LS(:,4), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5)
% histogram(scaled_WLS(:,4), 'Normalization', 'pdf','DisplayStyle','stairs', 'LineWidth',1.5, 'BinWidth', .5, 'LineStyle','--')
% title('Scaled Least Squares Errors')
% xlabel('Error (%)')
% ylabel('Probability Density')
% 
% legend('LS', 'WLS')

% Monte Carlo
% figure
% hold on
% plot(1:N, mean_LS_con,'LineWidth', 1.5)
% plot(1:N,mean_WLS_con,'LineWidth', 1.5)
% title('Mean convergence ')
% xlabel('Runs')
% ylabel('Mean (Absolute error)')
% legend('LS', 'WLS')
% 
% figure
% hold on
% plot(1:N,std_LS_con ,'LineWidth', 1.5)
% plot(1:N,std_WLS_con,'LineWidth', 1.5)
% title('std convergence ')
% xlabel('Runs')
% ylabel('STD (absolute error)')
% legend('LS', 'WLS')
% 
% Sorted Craters
% figure
% hold on
% plot(3:num_craters, mean_crat_LS(3:end)/1000,'LineWidth', 1.5)
% plot(3:num_craters,mean_crat_WLS(3:end)/1000,'LineWidth', 1.5)
% plot(3:num_craters,mean_crat_BCWLS(3:end)/1000,'g-','LineWidth', 1.5)
% title('Mean','FontSize',20,'FontName','Times New Roman')
% xlabel('Num Craters','FontSize',20,'FontName','Times New Roman')
% ylabel('Mean km (Absolute error)','FontSize',20,'FontName','Times New Roman')
% legend('LS', 'WLS','BCWLS','FontSize',20,'FontName','Times New Roman')
% grid on
% 
% figure
% hold on
% plot(3:num_craters,std_crat_LS(3:end)/1000 ,'LineWidth', 1.5)
% plot(3:num_craters,std_crat_WLS(3:end)/1000,'LineWidth', 1.5)
% plot(3:num_craters,std_crat_BCWLS(3:end)/1000,'g-','LineWidth', 1.5)
% title('std','FontSize',20,'FontName','Times New Roman')
% xlabel('Num Craters','FontSize',20,'FontName','Times New Roman')
% ylabel('STD km (absolute error)','FontSize',20,'FontName','Times New Roman')
% legend('LS', 'WLS','BCWLS','FontSize',20,'FontName','Times New Roman')
% grid on

altitudevec = altitude0:delta_altitude:altitudefinal;

% Variation over altitude/range

greylevels = linspace(0,0.8,5);

% mean
figure
hold on
for i = 1:length(solar_vec)
    shade = [greylevels(i) greylevels(i) greylevels(i)];
    plot(altitudevec/1000, mean_LS(i,:)/1000,'Color',shade,'LineStyle','-','LineWidth', 1.5)
    plot(altitudevec/1000,mean_WLS(i,:)/1000,'Color',shade,'LineStyle','--','LineWidth', 1.5)
    plot(altitudevec/1000,mean_BCWLS(i,:)/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
end
title('Mean','FontSize',20,'FontName','Times New Roman')
xlabel('Altitude km','FontSize',20,'FontName','Times New Roman')
ylabel('Mean km (Absolute error)','FontSize',20,'FontName','Times New Roman')
grid on
legend({'LS', 'WLS', 'BCWLS'},'FontSize',20,'FontName','Times New Roman')


figure % LS
hold on
for i = 1:length(solar_vec)
    shade = [greylevels(i) greylevels(i) greylevels(i)];
    plot(altitudevec/1000, mean_LS(i,:)/1000,'Color',shade,'LineStyle','-','LineWidth', 1.5)
end
title('Mean','FontSize',20,'FontName','Times New Roman')
xlabel('Altitude km','FontSize',20,'FontName','Times New Roman')
ylabel('Mean km (Absolute error)','FontSize',20,'FontName','Times New Roman')
grid on
legend({'LS'},'FontSize',20,'FontName','Times New Roman')

figure % WLS
hold on
for i = 1:length(solar_vec)
    shade = [greylevels(i) greylevels(i) greylevels(i)];
    plot(altitudevec/1000,mean_WLS(i,:)/1000,'Color',shade,'LineStyle','--','LineWidth', 1.5)
end
title('Mean','FontSize',20,'FontName','Times New Roman')
xlabel('Altitude km','FontSize',20,'FontName','Times New Roman')
ylabel('Mean km (Absolute error)','FontSize',20,'FontName','Times New Roman')
grid on
legend({'WLS'},'FontSize',20,'FontName','Times New Roman')

figure % BCWLS
hold on
for i = 1:length(solar_vec)
    shade = [greylevels(i) greylevels(i) greylevels(i)];
    plot(altitudevec/1000,mean_BCWLS(i,:)/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
end
title('Mean','FontSize',20,'FontName','Times New Roman')
xlabel('Altitude km','FontSize',20,'FontName','Times New Roman')
ylabel('Mean km (Absolute error)','FontSize',20,'FontName','Times New Roman')
grid on
legend({'BCWLS'},'FontSize',20,'FontName','Times New Roman')

% Standard Deviation
figure
hold on
for i = 1:length(solar_vec)
    shade = [greylevels(i) greylevels(i) greylevels(i)];
    plot(altitudevec/1000,std_LS(i,:) /1000,'Color',shade,'LineStyle','-','LineWidth', 1.5)
    plot(altitudevec/1000,std_WLS(i,:)/1000,'Color',shade,'LineStyle','--','LineWidth', 1.5)
    plot(altitudevec/1000,std_BCWLS(i,:)/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
end
title('std','FontSize',20,'FontName','Times New Roman')
xlabel('Altitude km','FontSize',20,'FontName','Times New Roman')
ylabel('STD km (absolute error)','FontSize',20,'FontName','Times New Roman')
grid on
legend({'WLS'},'FontSize',20,'FontName','Times New Roman')