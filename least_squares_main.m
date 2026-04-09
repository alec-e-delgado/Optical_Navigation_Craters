close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEAST SQUARES (LS) SCRIPT
% This script computes various LS estimates for the Spacecraft (SC)
% position and compares their absolute errors
%
% DESCRIPTION
%  - Defines simulation constants and initial conditions
%  - Calls dynamics, control, and plotting functions
%
% EXPERIMENTS
%  - Convergence: Tests for the number of Monte Carlo Simulations needed
%  for the mean error and std to converge
%  - Altitudes: Estimates poisitional error over various altitudes
%
% INPUTS
%  - Which experiments to run
%  - Experiment inputs
%
% OUTPUTS
%  - Plotted simulation results
%
% Author: Alec Delgado
% Date: 2026-03-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT

% Define which tests to run (true or false)
test_convergence = true;
test_altitudes = true;

% Convergence Experiment
distance_conv = 35000;    % (km) run at this altitude
N = 10000;                % runs for Monte Carlo Simulation

% Altitudes Experiment
ditance0 = 35000;         % (km) initial altitude above moons surface 
distance_final = 97000;   % (km) final altitutde to estimate at  
delta_distance = 500;     % (km) 
iter = 1;                 % iterator


%% Constants 
moon_angle = 90;                        % (deg) incidence angle with moon surface 
radius_Moon = 1.7374e6;                 % (m)
inertial_uv = [0.4811; 0.2580; 0.8379]; % direction of spacecraft position

%% Convergence Experiment

if test_convergence % run if desired

    fprintf("------------------------ Convergence Test ------------------------\n")

    distance_sc = distance_conv*10^3;                % (m) Initialize 
    r_sc_I = inertial_uv*(distance_sc); % (m) Inertial position vector

    % Get detected crater uncertainty
    [~,repeat_matrix_detections,~,~,~] = angular_error_calc(distance_sc/1000, moon_angle);
    mat_size = size(repeat_matrix_detections);
    num_craters = mat_size(1); 

    fprintf('Detected %d craters\n\n', num_craters);

    % Extract statistical quantities from matrix for each crater
    sc_bearings = repeat_matrix_detections(:, 3);   % (deg) angle from normal to crater
    std_dev_norm = repeat_matrix_detections(:,2);   % (normalized) standard deviation
    mean_norm = repeat_matrix_detections(:,1);      % (normalized) mean
    crater_radius = repeat_matrix_detections(:, 4); % (m)

    % Un-normalize standard deviation and mean values by moon angular area
    moon_ang_area = 2*rad2deg(asin(1737.4/(distance_sc/1000))); % normalization factor
    std_dev_unnorm = deg2rad(std_dev_norm*moon_ang_area/100);   % unnormalized
    mean_unnorm = deg2rad(mean_norm*moon_ang_area/100);         % unnormalized

    % Treat ill standard deviation values
    idx_Nan = find(isnan(std_dev_unnorm)); % finds indices of Nan values
    std_dev_unnorm(idx_Nan) = [];          % remove Nan values
    mean_unnorm(idx_Nan) = [];
    sc_bearings(idx_Nan) = [];
    crater_radius(idx_Nan) = [];

    idx_zero = find(std_dev_unnorm == 0);  % finds indices of 0 std devs value
    std_dev_unnorm(idx_zero) = [];         % remove 0 values
    mean_unnorm(idx_zero) = [];
    sc_bearings(idx_zero) = [];
    crater_radius(idx_zero) = [];
    
    % Initialize errors 
    conv_LS_errors = zeros(1, N);
    conv_mean_LS = zeros(1, N);
    conv_std_LS = zeros(1, N);

    conv_WLS_errors = zeros(1, N);
    conv_mean_WLS = zeros(1, N);
    conv_std_WLS = zeros(1, N);

    conv_BCWLS_errors = zeros(1, N);
    conv_mean_BCWLS = zeros(1, N);
    conv_std_BCWLS = zeros(1, N);

    for i = 1:N % Monte Carlo iterations
        % Generate angular measurements and get uncertainty
        [std_elevation,std_azimuth,meas_azimuth,meas_elevation,r_crater_LOS_nom] = crater_pos(r_sc_I,sc_bearings,std_dev_unnorm);

        % LS estimate error
        [pos_LS_est,LS_error,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation,r_crater_LOS_nom); 

        % WLS estimate error
        [pos_WLS_est,WLS_error,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,r_crater_LOS_nom,pos_LS_est,A_mat,z_mat);

        % BCWLS esimtate error
        [pos_BCWLS_est,BCWLS_error] = BCWLS_estimate(std_azimuth,std_elevation,meas_azimuth,meas_elevation,pos_WLS_est,A_mat,W_block,do);


        % Store values
        conv_LS_errors(i) = LS_error;
        conv_WLS_errors(i) = WLS_error;
        conv_BCWLS_errors(i) = BCWLS_error;

        conv_mean_LS(i) = mean(conv_LS_errors(1:i));
        conv_std_LS(i) = std(conv_LS_errors(1:i),0,2);

        conv_mean_WLS(i) = mean(conv_WLS_errors(1:i));
        conv_std_WLS(i) = std(conv_WLS_errors(1:i),0,2);

        conv_mean_BCWLS(i) = mean(conv_BCWLS_errors(1:i));
        conv_std_BCWLS(i) = std(conv_BCWLS_errors(1:i),0,2);
    
    end

    % Plots

    figure % mean
    hold on
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Times New Roman';
    plot(1:N,conv_mean_LS/1000,'r-','LineWidth', 1.5)
    plot(1:N,conv_mean_WLS/1000,'k-','LineWidth', 1.5)
    plot(1:N,conv_mean_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Runs')
    ylabel('Mean Range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',18)
    
    figure % standard deviation
    hold on
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Times New Roman';
    plot(1:N,conv_std_LS/1000 ,'r-','LineWidth', 1.5)
    plot(1:N,conv_std_WLS/1000,'k-','LineWidth', 1.5)
    plot(1:N,conv_std_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Runs')
    ylabel('STD Range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',18)
end

%% Altitudes Experiment 

if test_altitudes

    fprintf("------------------------- Altitudes Test -------------------------\n")

    distance_sc = ditance0*10^3;        % (m) initialize
    r_sc_I = inertial_uv*(distance_sc); % (m) Inertial position vector

    num_iter = round((distance_final-ditance0)/delta_distance) + 1; % for initialization
    dist_list = zeros(num_iter,1);                                  % storing distances

    % Initialize errors 
    dist_LS_errors = zeros(num_iter, N); % [num_iterxN]
    dist_mean_LS = zeros(num_iter, 1);   % [num_iterx1]
    dist_std_LS = zeros(num_iter, 1);

    dist_WLS_errors = zeros(num_iter, N);
    dist_mean_WLS = zeros(num_iter, 1);
    dist_std_WLS = zeros(num_iter, 1);

    dist_BCWLS_errors = zeros(num_iter, N);
    dist_mean_BCWLS = zeros(num_iter, 1);
    dist_std_BCWLS = zeros(num_iter, 1);

    while distance_sc <= distance_final*10^3 % m

        fprintf("Computing estimate at %.f km\n", distance_sc/1000)

        dist_list(iter) = distance_sc; % m

        % Get detected crater uncertainty
        [~,repeat_matrix_detections,~,~,~] = angular_error_calc(distance_sc/1000, moon_angle);
        mat_size = size(repeat_matrix_detections);
        num_craters = mat_size(1); 
    
        % Extract statistical quantities from matrix for each crater
        sc_bearings = repeat_matrix_detections(:, 3);   % (deg) angle from normal to crater
        std_dev_norm = repeat_matrix_detections(:,2);   % (normalized) standard deviation
        mean_norm = repeat_matrix_detections(:,1);      % (normalized) mean
        crater_radius = repeat_matrix_detections(:, 4); % (m)
    
        % Un-normalize standard deviation and mean values by moon angular area
        moon_ang_area = 2*rad2deg(asin(1737.4/(distance_sc/1000))); % normalization factor
        std_dev_unnorm = deg2rad(std_dev_norm*moon_ang_area/100);   % unnormalized
        mean_unnorm = deg2rad(mean_norm*moon_ang_area/100);         % unnormalized
    
        % Treat ill standard deviation values
        idx_Nan = find(isnan(std_dev_unnorm)); % finds indices of Nan values
        std_dev_unnorm(idx_Nan) = [];          % remove Nan values
        mean_unnorm(idx_Nan) = [];
        sc_bearings(idx_Nan) = [];
        crater_radius(idx_Nan) = [];
    
        idx_zero = find(std_dev_unnorm == 0);  % finds indices of 0 std devs value
        std_dev_unnorm(idx_zero) = [];         % remove 0 values
        mean_unnorm(idx_zero) = [];
        sc_bearings(idx_zero) = [];
        crater_radius(idx_zero) = [];

        for i = 1:N % Monte Carlo iterations

            % Generate angular measurements and get uncertainty
            [std_elevation,std_azimuth,meas_azimuth,meas_elevation,r_crater_LOS_nom] = crater_pos(r_sc_I,sc_bearings,std_dev_unnorm);
    
            % LS estimate error
            [pos_LS_est,LS_error,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation,r_crater_LOS_nom); 
    
            % WLS estimate error
            [pos_WLS_est,WLS_error,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,r_crater_LOS_nom,pos_LS_est,A_mat,z_mat);
    
            % BCWLS esimtate error
            [pos_BCWLS_est,BCWLS_error] = BCWLS_estimate(std_azimuth,std_elevation,meas_azimuth,meas_elevation,pos_WLS_est,A_mat,W_block,do);
    
            % Store values
            dist_LS_errors(iter,i) = LS_error;
            dist_WLS_errors(iter,i) = WLS_error;
            dist_BCWLS_errors(iter,i) = BCWLS_error;
        end

        % Calculate Mean and STD values
        dist_mean_LS(iter) = mean(dist_LS_errors(iter,:));         % LS
        dist_std_LS(iter) = std(dist_LS_errors(iter,:),0,2);

        dist_mean_WLS(iter) = mean(dist_WLS_errors(iter,:));       % WLS
        dist_std_WLS(iter) = std(dist_WLS_errors(iter,:),0,2);

        dist_mean_BCWLS(iter) = mean(dist_BCWLS_errors(iter,:));   % BCWLS
        dist_std_BCWLS(iter) = std(dist_BCWLS_errors(iter,:),0,2);

        % Update Values
        distance_sc = distance_sc + delta_distance*10^3; % (m)
        iter = iter + 1;
    end

    % Plots

    figure % mean
    hold on
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Times New Roman';
    plot(dist_list/1000000,dist_mean_LS/1000,'r-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_mean_WLS/1000,'k-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_mean_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Distance (km $\times$ 10$^3$) ','Interpreter','latex')
    ylabel('Mean Range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',18)
    
    figure % standard deviation
    hold on
    ax = gca;
    ax.FontSize = 18;
    ax.FontName = 'Times New Roman';
    plot(dist_list/1000000,dist_std_LS/1000 ,'r-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_std_WLS/1000,'k-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_std_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Distance (km $\times$ 10$^3$)','Interpreter','latex')
    ylabel('STD Range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',18)
end



% while altitude_sc <= altitudefinal

% Number of times to run simulation

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
conv_LS_errors = zeros(1, N);
conv_std_LS = zeros(1, N);
conv_mean_WLS = zeros(1, N);
std_WLS_con = zeros(1, N);


% for i = 1:N % Monte Carlo
% for i = 1 % need minimum three measurements
    for i = 1:N % num iterations
    % Pass statistic values 
    [BCWLS_e, WLS_e, LS_error, r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearings, std_dev_unnorm, mean_unnorm,crater_radius);
    
        LS_errors(:,i) = LS_error(end); % end only for altitude variation
        WLS_errors(:,i) = WLS_e(end);
        BCWLS_errors(:,i) = BCWLS_e(end);

        crater_LS(i,:) = LS_error;
        crater_WLS(i,:) = WLS_e;
        crater_BCWLS(i,:) = BCWLS_e;

        % mean_LS_con(j) = mean(LS_errors(end,1:j),2);
        % std_LS_con(j) = std(LS_errors(end,1:j),0,2);
        % 
        % mean_WLS_con(j) = mean(WLS_errors(end,1:j),2);
        % std_WLS_con(j) = std(WLS_errors(end,1:j),0,2);
    
    end
mean_LS(z) = mean(LS_errors,2);
std_LS(z) = std(LS_errors,0,2);
mean_crat_LS = mean(crater_LS,1);
std_crat_LS = std(crater_LS,0,1);

mean_WLS(z) = mean(WLS_errors,2);
std_WLS(z) = std(WLS_errors,0,2);
mean_crat_WLS = mean(crater_WLS,1);
std_crat_WLS = std(crater_WLS,0,1);

mean_BCWLS(z) = mean(BCWLS_errors,2);
std_BCWLS(z) = std(BCWLS_errors,0,2);
mean_crat_BCWLS = mean(crater_BCWLS,1);
std_crat_BCWLS = std(crater_BCWLS,0,1);

fprintf("Iteration: %.d\n", distance_sc)

% Update iterators
distance_sc = distance_sc + delta_distance; % add 500 m each iteration
z = z + 1;
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
figure
hold on
plot(3:num_craters, mean_crat_LS(3:end)/1000,'LineWidth', 1.5)
plot(3:num_craters,mean_crat_WLS(3:end)/1000,'LineWidth', 1.5)
plot(3:num_craters,mean_crat_BCWLS(3:end)/1000,'g-','LineWidth', 1.5)
title('Mean','FontSize',20,'FontName','Times New Roman')
xlabel('Num Craters','FontSize',20,'FontName','Times New Roman')
ylabel('Mean km (Absolute error)','FontSize',20,'FontName','Times New Roman')
legend('LS', 'WLS','BCWLS','FontSize',20,'FontName','Times New Roman')
grid on

figure
hold on
plot(3:num_craters,std_crat_LS(3:end)/1000 ,'LineWidth', 1.5)
plot(3:num_craters,std_crat_WLS(3:end)/1000,'LineWidth', 1.5)
plot(3:num_craters,std_crat_BCWLS(3:end)/1000,'g-','LineWidth', 1.5)
title('std','FontSize',20,'FontName','Times New Roman')
xlabel('Num Craters','FontSize',20,'FontName','Times New Roman')
ylabel('STD km (absolute error)','FontSize',20,'FontName','Times New Roman')
legend('LS', 'WLS','BCWLS','FontSize',20,'FontName','Times New Roman')
grid on

altitudevec = ditance0:delta_distance:distance_final;

% Variation over altitude/range
% figure
% hold on
% plot(altitudevec, mean_LS,'LineWidth', 1.5)
% plot(altitudevec,mean_WLS,'LineWidth', 1.5)
% plot(altitudevec,mean_BCWLS,'LineWidth',1.5)
% title('Mean')
% xlabel('Altitude')
% ylabel('Mean (Absolute error)')
% legend('LS', 'WLS','BCWLS')
% 
% figure
% hold on
% plot(altitudevec,std_LS ,'LineWidth', 1.5)
% plot(altitudevec,std_WLS,'LineWidth', 1.5)
% plot(altitudevec,std_BCWLS,'LineWidth',1.5)
% title('std')
% xlabel('Altitude')
% ylabel('STD (absolute error)')
% legend('LS', 'WLS','BCWLS')