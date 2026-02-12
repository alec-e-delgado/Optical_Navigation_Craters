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

% Number of times to run simulation

N = 10000; % to use for Monte Carlo

% Pre allocate space for errors
mean_LS = zeros(N, 1);
std_LS = zeros(N, 1);

mean_WLS = zeros(N, 1);
std_WLS = zeros(N, 1);

% for j = 1:jend

leave = 1;

[azimuth, elevation, r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearings, std_devs, means,leave);
[azimuth_sort, idx_sort] = sort(azimuth,'ascend');
sc_bearings = sc_bearings(idx_sort);
std_devs = std_devs(idx_sort);
means = means(idx_sort);

% Need to put above implementation in crater_pos, organize and just
% extracxt i number of craters

% The azimuth and elevation angles are not the same for each call of the
% function, thats why it needs to be added within the function, just add it
% there, sort, and extract the best i craters, use i as input to function

% Pre allocate space
LS_errors = zeros(num_craters, N); % Preallocate LS_errors matrix
WLS_errors = zeros(num_craters, N); % Preallocate WLS_errors matrix
mean_LS_con = zeros(1, N);
std_LS_con = zeros(1, N);
mean_WLS_con = zeros(1, N);
std_WLS_con = zeros(1, N);


% for i = 1:N % Monte Carlo
% for i = 1 % need minimum three measurements
    for j = 1:N % num iterations
    % Pass statistic values 
    [WLS_e, LS_e, r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearings, std_devs, means,leave);
    
        LS_errors(:,j) = LS_e;
        WLS_errors(:,j) = WLS_e;
        fprintf("Iteration: %.d\n", j)

        mean_LS_con(j) = mean(LS_errors(end,1:j),2);
        std_LS_con(j) = std(LS_errors(end,1:j),0,2);

        mean_WLS_con(j) = mean(WLS_errors(end,1:j),2);
        std_WLS_con(j) = std(WLS_errors(end,1:j),0,2);
    
    end
mean_LS = mean(LS_errors,2);
std_LS = std(LS_errors,0,2);

mean_WLS = mean(WLS_errors,2);
std_WLS = std(WLS_errors,0,2);
% end


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

% Monte Carlo
figure
hold on
plot(1:N, mean_LS_con,'LineWidth', 1.5)
plot(1:N,mean_WLS_con,'LineWidth', 1.5)
title('Mean convergence ')
xlabel('Runs')
ylabel('Mean (Absolute error)')
legend('LS', 'WLS')

figure
hold on
plot(1:N,std_LS_con ,'LineWidth', 1.5)
plot(1:N,std_WLS_con,'LineWidth', 1.5)
title('std convergence ')
xlabel('Runs')
ylabel('STD (absolute error)')
legend('LS', 'WLS')

% Sorted Craters
figure
hold on
plot(3:num_craters, mean_LS(3:end),'LineWidth', 1.5)
plot(3:num_craters,mean_WLS(3:end),'LineWidth', 1.5)
title('Mean')
xlabel('Num Craters')
ylabel('Mean (Absolute error)')
legend('LS', 'WLS')

figure
hold on
plot(3:num_craters,std_LS(3:end) ,'LineWidth', 1.5)
plot(3:num_craters,std_WLS(3:end),'LineWidth', 1.5)
title('std')
xlabel('Num Craters')
ylabel('STD (absolute error)')
legend('LS', 'WLS')
