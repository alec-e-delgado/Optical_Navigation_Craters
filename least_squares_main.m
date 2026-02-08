close all; clear; clc;

% The least squares solution will take n number of measurements (at least
% 3) to compute the position of the spacecraft

% Define parameters
altitude_sc = 35000; % current distance away from moon surface (need units)
moon_angle = 90; % incidence angle with moon surface (degrees)
radius_Moon = 1.7374e6; % m
distance_sc = altitude_sc+radius_Moon; % m

% Calculate position vector of spacecraft used to locate craters using some
% random unit vector
r_sc_I = [0.4811; 0.2580; 0.8379]*(distance_sc); % m

% Calculate matrix of mean and std
[angular_errors,repeat_matrix_detections,x_2, y_2, z_2] = angular_error_calc(altitude_sc, moon_angle);

% Max index value
rep_matrix_size = size(repeat_matrix_detections);
fprintf('Max index: %d\n', rep_matrix_size(1));
is_pass = false; % used to ensure index is chosen within the size of the array

% Display list of craters
% for j = 1:rep_matrix_size(1)
%     fprintf('%d -- sc_inc = %.2f -- crater_radius = %.2f\n', j, repeat_matrix_detections(j,3),repeat_matrix_detections(j,4))
% end

desired_indices = 1:length(x_2); % list of indices to use for crater generation

% Extract crater quantities associated with desired indeces
std_devs = repeat_matrix_detections(desired_indices,2); % rad
means = repeat_matrix_detections(desired_indices,1); % rad
sc_bearings = repeat_matrix_detections(desired_indices, 3); % rad

% Number of times to run simulation

N = 10000;

% Pre allocate space for errors
LS_errors = zeros(N, 1);
WLS_errors = zeros(N, 1);

mean_LS = zeros(N, 1);
std_LS = zeros(N, 1);

mean_WLS = zeros(N, 1);
std_WLS = zeros(N, 1);

% for j = 1:jend

for i = 1:N
% Pass statistic values 
    [WLS_e, LS_e, r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearings, std_devs, means);

    LS_errors(i) = LS_e;
    WLS_errors(i) = WLS_e;

    mean_LS(i) = mean(LS_errors(1:i));
    std_LS(i) = std(LS_errors(1:i));

    mean_WLS(i) = mean(WLS_errors(1:i));
    std_WLS(i) = std(WLS_errors(1:i));
end

% % Scale errors by altitude
% scaled_LS(:,j) = LS_errors(:,j)/altitude_sc * 100; 
% scaled_WLS(:,j) = WLS_errors(:,j)/altitude_sc * 100;
% end



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

figure
hold on
plot(1:N, mean_LS,'LineWidth', 1.5)
plot(1:N,mean_WLS,'LineWidth', 1.5)
title('Mean convergence ')
xlabel('Runs')
ylabel('Mean (Absolute error)')
legend('LS', 'WLS')

figure
hold on
plot(1:N,std_LS ,'LineWidth', 1.5)
plot(1:N,std_WLS,'LineWidth', 1.5)
title('std convergence ')
xlabel('Runs')
ylabel('STD (absolute error)')
legend('LS', 'WLS')



