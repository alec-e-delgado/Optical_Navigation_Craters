close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEAST SQUARES (LS) SCRIPT
% This script computes various LS estimates for the Spacecraft (SC)
% position and compares their absolute errors
%
% DESCRIPTION
%  - Runs several experiments to characterize LS, WLS, BCWLS estimators 
%  from angle only measurements
%
% EXPERIMENTS
%  - Angular error: Shows histogram of angular errors for bearing, azimuth,
%    and elevation angles
%  - Convergence: Tests for the number of Monte Carlo Simulations needed
%    for the mean error and std to converge
%       - Includes histogram plots
%  - Distances: Estimates poisitional error over various distances
%  - Solar Phase: Altitude test across different solar phase angles, only
%    considers BCWLS
%       - Computed using BCWLS estimate
%  - Sorting: Sort the measurements using different metrics and use the n
%    best
%       - angular std, crater radius, azimuth
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
test_angular_error = true;
test_convergence = false;
test_distances = false;
test_solar_phase = false;
test_sort = false; 

% Convergence Experiment
distance_conv = 10000;    % (km) run at this altitude
N = 10000;                % runs for Monte Carlo Simulation

% Distances Experiment
ditance0 = 10000;         % (km) initial distance from Moon center
distance_final = 97000;   % (km) final distance to estimate at  
delta_distance = 1000;    % (km) 
iter = 1;                 % iterator

% Solar Phase Experiment
solar_angle0 = 0;         % (deg) initial solar angle
solar_angle_final = 120;  % (deg) final solar angle 
solar_delta = 20;         % (deg) 

% Sorting Experiment
distance_sorting = 45000; % (km)

%% Constants 
moon_angle = 90;                        % (deg) incidence angle with moon surface 
radius_Moon = 1.7374e6;                 % (m)
inertial_uv = [0.4811; 0.2580; 0.8379]; % direction of spacecraft position

%% Angular Error Experiment

if test_angular_error

    fprintf("------------------------ Angular Error Test ------------------------\n")

    distance_sc = distance_conv*10^3;   % (m) Initialize 
    r_sc_I = inertial_uv*(distance_sc); % (m) Inertial position vector

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
    std_dev_unnorm = deg2rad(std_dev_norm*moon_ang_area/100);   % (rad) unnormalized
    mean_unnorm = deg2rad(mean_norm*moon_ang_area/100);         % (rad) unnormalized

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

    % Select one measurement and expand for visualization
    bearing_test = sc_bearings(1);
    std_test = std_dev_unnorm(1);

    % Initialize
    bearing_error_vec = zeros(1,N);
    azimuth_error_vec = zeros(1,N);
    elevation_error_vec = zeros(1,N);
    
    for i = 1:N
        [bearing_error,azimuth_error,elevation_error,std_elevation,std_azimuth] = test(r_sc_I,bearing_test,std_test);
        bearing_error_vec(i) = bearing_error;
        azimuth_error_vec(i) = azimuth_error;
        elevation_error_vec(i) = elevation_error;
    end

    % Print standard deviations of all angles
    fprintf('STD Bearing angle:   %.5f deg\n', rad2deg(std_test))     
    fprintf('STD Azimuth angle:   %.5f deg\n', rad2deg(std_azimuth))
    fprintf('STD Elevation angle: %.5f deg\n', rad2deg(std_elevation))
    
    figure
    histogram(rad2deg(bearing_error_vec),'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5)
    box on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    xlabel('Bearing angle measurement error (deg)')
    ylabel('Probability density (%)')

    figure(100)
    subplot(1,2,1); histogram(rad2deg(azimuth_error_vec),'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5)
    box on; xlabel('Azmiuth measurement error (deg)','FontName','Times New Roman','FontSize',14)
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    subplot(1,2,2); histogram(rad2deg(elevation_error_vec),'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5)
    box on; xlabel('Elevation measurement error (deg)','FontName','Times New Roman','FontSize',14)
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    han=axes(figure(100),'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Probability density (%)');
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';


end

%% Convergence Experiment

if test_convergence % run if desired

    fprintf("------------------------ Convergence Test ------------------------\n")

    distance_sc = distance_conv*10^3;   % (m) Initialize 
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
    std_dev_unnorm = deg2rad(std_dev_norm*moon_ang_area/100);   % (rad) unnormalized
    mean_unnorm = deg2rad(mean_norm*moon_ang_area/100);         % (rad) unnormalized

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

    range_error_LS = zeros(1,N); % for histogram generation
    range_error_WLS = zeros(1,N);
    range_error_BCWLS = zeros(1,N);

    for i = 1:N % Monte Carlo iterations
        % Generate angular measurements and get uncertainty
        [std_elevation,std_azimuth,meas_azimuth,meas_elevation,r_crater_LOS_nom,R_LOS_I] = crater_pos(r_sc_I,sc_bearings,std_dev_unnorm);

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

        % Rotate position estimates for histogram generation
        pos_inertial_LS = R_LOS_I*(pos_LS_est-[0;distance_sc;0]);
        pos_inertial_WLS = R_LOS_I*(pos_WLS_est-[0;distance_sc;0]);
        pos_inertial_BCWLS = R_LOS_I*(pos_BCWLS_est-[0;distance_sc;0]);

        range_error_LS(i) = norm(pos_inertial_LS) - norm(r_sc_I);
        range_error_WLS(i) = norm(pos_inertial_WLS) - norm(r_sc_I);
        range_error_BCWLS(i) = norm(pos_inertial_BCWLS) - norm(r_sc_I);
    
    end
    
    % Calculatre standard deviations of each estimate method
    std_LS = std(range_error_LS/1000);       % (km)
    std_WLS = std(range_error_WLS/1000);
    std_BCWLS = std(range_error_BCWLS/1000); 
    
    % Print standard deviations
    fprintf('LS STD:    %.3f km \n',std_LS);
    fprintf('WLS STD:   %.3f km \n',std_WLS);
    fprintf('BCWLS STD: %.3f km \n',std_BCWLS);

    % Plots

    figure % mean
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(1:N,conv_mean_LS/1000,'r-','LineWidth', 1.5)
    plot(1:N,conv_mean_WLS/1000,'k-','LineWidth', 1.5)
    plot(1:N,conv_mean_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Runs')
    ylabel('Mean absolute range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',14)
    box on
    
    figure % standard deviation
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(1:N,conv_std_LS/1000 ,'r-','LineWidth', 1.5)
    plot(1:N,conv_std_WLS/1000,'k-','LineWidth', 1.5)
    plot(1:N,conv_std_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Runs')
    ylabel('STD absolute range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',14)
    box on

    figure % range error
    hold on
    histogram(range_error_LS/1000,'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5,'BinWidth',.200,'EdgeColor', 'r')
    histogram(range_error_WLS/1000,'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5,'BinWidth',.200,'EdgeColor', 'k')
    histogram(range_error_BCWLS/1000,'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5,'BinWidth',.200,'EdgeColor', 'g','LineStyle','--')
    box on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    xlabel('Range Error (km)')
    ylabel('Probability density (%)')
    legend('LS','WLS','BCWLS','FontName', 'Times New Roman', 'FontSize', 14)

    figure % absolute range error
    hold on
    histogram(conv_LS_errors/1000,'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5,'BinWidth',.200,'EdgeColor', 'r')
    histogram(conv_WLS_errors/1000,'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5,'BinWidth',.200,'EdgeColor', 'k')
    histogram(conv_BCWLS_errors/1000,'Normalization','percentage','DisplayStyle','stairs', 'LineWidth',1.5,'BinWidth',.200,'EdgeColor', 'g','LineStyle','--')
    box on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    xlabel('Absolute Range Error (km)')
    ylabel('Probability density (%)')
    legend('LS','WLS','BCWLS','FontName', 'Times New Roman', 'FontSize', 14)
end

%% Distances Experiment 

if test_distances

    fprintf("------------------------- Distances Test -------------------------\n")
    
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

    % Coordiante estimate errors in LOS reference frame
    pos_error_LOS = zeros(num_iter, N); % (y) Line-of-sight axis 
    pos_error_HCA = zeros(num_iter, N); % (x) Horizontal cross-axis
    pos_error_VCA = zeros(num_iter, N); % (z) Vertical cross-axis

    pos_mean_LOS = zeros(num_iter, 1);  % (m) average coordinate errors
    pos_mean_HCA = zeros(num_iter, 1); 
    pos_mean_VCA = zeros(num_iter, 1); 

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
            [std_elevation,std_azimuth,meas_azimuth,meas_elevation,r_crater_LOS_nom,R_LOS_I] = crater_pos(r_sc_I,sc_bearings,std_dev_unnorm);
    
            % LS estimate error
            [pos_LS_est,LS_error,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation,r_crater_LOS_nom); 
    
            % WLS estimate error
            [pos_WLS_est,WLS_error,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,r_crater_LOS_nom,pos_LS_est,A_mat,z_mat);
    
            % BCWLS esimtate error
            [pos_BCWLS_est,BCWLS_error] = BCWLS_estimate(std_azimuth,std_elevation,meas_azimuth,meas_elevation,pos_WLS_est,A_mat,W_block,do);
    
            % Store values
            dist_LS_errors(iter,i) = LS_error;               % (m) Range errors
            dist_WLS_errors(iter,i) = WLS_error;
            dist_BCWLS_errors(iter,i) = BCWLS_error;

            pos_error_LOS(iter,i) = abs(pos_BCWLS_est(2));   % (m) Coordinate errors
            pos_error_HCA(iter,i) = abs(pos_BCWLS_est(1));
            pos_error_VCA(iter,i) = abs(pos_BCWLS_est(3));
            
        end

        % Calculate Mean and STD values
        dist_mean_LS(iter) = mean(dist_LS_errors(iter,:));         % LS
        dist_std_LS(iter) = std(dist_LS_errors(iter,:),0,2);

        dist_mean_WLS(iter) = mean(dist_WLS_errors(iter,:));       % WLS
        dist_std_WLS(iter) = std(dist_WLS_errors(iter,:),0,2);

        dist_mean_BCWLS(iter) = mean(dist_BCWLS_errors(iter,:));   % BCWLS
        dist_std_BCWLS(iter) = std(dist_BCWLS_errors(iter,:),0,2);

        pos_mean_LOS(iter) = mean(pos_error_LOS(iter,:)); 
        pos_mean_HCA(iter) = mean(pos_error_HCA(iter,:));
        pos_mean_VCA(iter) = mean(pos_error_VCA(iter,:));

        % Update Values
        distance_sc = distance_sc + delta_distance*10^3; % (m)
        r_sc_I = inertial_uv*(distance_sc); % (m) Inertial position vector
        iter = iter + 1;
    end

    % Plots

    figure % mean
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(dist_list/1000000,dist_mean_LS/1000,'r-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_mean_WLS/1000,'k-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_mean_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Distance (km $\times$ 10$^3$) ','Interpreter','latex')
    ylabel('Mean absolute range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',14)
    box on;
    
    figure % standard deviation
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(dist_list/1000000,dist_std_LS/1000 ,'r-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_std_WLS/1000,'k-','LineWidth', 1.5)
    plot(dist_list/1000000,dist_std_BCWLS/1000,'g--','LineWidth', 1.5)
    xlabel('Distance (km $\times$ 10$^3$)','Interpreter','latex')
    ylabel('STD absolute range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',14)
    box on;

    fig = figure; % Coordinate errors
    subplot(3,1,1); plot(dist_list/1000000,pos_mean_LOS/1000,'g-','LineWidth', 1.5)
    ax = gca; ax.FontSize = 14; ax.FontName = 'Times New Roman';
    subplot(3,1,2); plot(dist_list/1000000,pos_mean_HCA/1000,'g-','LineWidth', 1.5)
    ax = gca; ax.FontSize = 14; ax.FontName = 'Times New Roman';
    subplot(3,1,3); plot(dist_list/1000000,pos_mean_VCA/1000,'g-','LineWidth', 1.5)
    ax = gca; ax.FontSize = 14; ax.FontName = 'Times New Roman';
    han=axes(fig,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Absolute position error (km)');
    xlabel(han,'Distance (km $\times$ 10$^3$)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    box on;

    figure % difference between WLS and BCWLS
    plot(dist_list/1000000,dist_mean_LS/1000,'LineWidth', 1.5)
    xlabel('Distance (km $\times$ 10$^3$) ','Interpreter','latex')
    ylabel('$\mu_{WLS}-\mu_{BCWLS}$ (m)','Interpreter','latex')
    ax = gca; ax.FontSize = 14; ax.FontName = 'Times New Roman';
    box on;

end

%% Solar Phase Experiment

if test_solar_phase

    fprintf("------------------------- Solar Phase Test -------------------------\n")
    
    num_solar_phase = (solar_angle_final-solar_angle0)/solar_delta + 1;
    solar_angle = solar_angle0;                   % (deg) Initialize
    solar_idx = 1;                                % index for plotting

    % Create shades of green for plotting (dark -> light)
    greens = [0 71 0; % dark
              0 117 0;
              0 163 0;
              0 209 0;
              0 255 0;
              138 255 138;
              184 255 184]/255; % light
    
    while solar_angle <= solar_angle_final 

    fprintf("Computing estimate at %.f deg\n", solar_angle)
    
    distance_sc = ditance0*10^3;        % (m) initialize
    r_sc_I = inertial_uv*(distance_sc); % (m) Inertial position vector
    iter = 1;

    num_iter = round((distance_final-ditance0)/delta_distance) + 1; % for initialization
    dist_list = zeros(num_iter,1);                                  % storing distances

    % Initialize errors 
    dist_BCWLS_errors = zeros(num_iter, N);
    dist_mean_BCWLS = zeros(num_iter, 1);
    dist_std_BCWLS = zeros(num_iter, 1);


    % Coordiante estimate errors in LOS reference frame
    pos_error_LOS = zeros(num_iter, N); % (y) Line-of-sight axis 
    pos_error_HCA = zeros(num_iter, N); % (x) Horizontal cross-axis
    pos_error_VCA = zeros(num_iter, N); % (z) Vertical cross-axis

    pos_mean_LOS = zeros(num_iter, 1);  % (m) average coordinate errors
    pos_mean_HCA = zeros(num_iter, 1); 
    pos_mean_VCA = zeros(num_iter, 1);

    while distance_sc <= distance_final*10^3 % m

        dist_list(iter) = distance_sc; % m

        % Get detected crater uncertainty
        [~,repeat_matrix_detections,~,~,~] = angular_error_calc(distance_sc/1000, solar_angle);
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
            [std_elevation,std_azimuth,meas_azimuth,meas_elevation,r_crater_LOS_nom,R_LOS_I] = crater_pos(r_sc_I,sc_bearings,std_dev_unnorm);
    
            % LS estimate error
            [pos_LS_est,LS_error,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation,r_crater_LOS_nom); 
    
            % WLS estimate error
            [pos_WLS_est,WLS_error,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,r_crater_LOS_nom,pos_LS_est,A_mat,z_mat);
    
            % BCWLS esimtate error
            [pos_BCWLS_est,BCWLS_error] = BCWLS_estimate(std_azimuth,std_elevation,meas_azimuth,meas_elevation,pos_WLS_est,A_mat,W_block,do);
    
            % Store values
            dist_BCWLS_errors(iter,i) = BCWLS_error; % (m) range error

            pos_error_LOS(iter,i) = abs(pos_BCWLS_est(2));   % (m) Coordinate errors
            pos_error_HCA(iter,i) = abs(pos_BCWLS_est(1));
            pos_error_VCA(iter,i) = abs(pos_BCWLS_est(3));
        end

        % Calculate Mean and STD values
        dist_mean_BCWLS(iter) = mean(dist_BCWLS_errors(iter,:));   % BCWLS
        dist_std_BCWLS(iter) = std(dist_BCWLS_errors(iter,:),0,2);

        pos_mean_LOS(iter) = mean(pos_error_LOS(iter,:));          % position errors
        pos_mean_HCA(iter) = mean(pos_error_HCA(iter,:));
        pos_mean_VCA(iter) = mean(pos_error_VCA(iter,:));

        % Update Values
        distance_sc = distance_sc + delta_distance*10^3; % (m)
        r_sc_I = inertial_uv*(distance_sc);              % (m) Inertial position vector
        iter = iter + 1;
    end

    % Update Solar Angle 
    solar_angle = solar_angle + solar_delta; % (deg)
    
    % Plots
    shade = [greens(solar_idx,1) greens(solar_idx,2) greens(solar_idx,3)];
    
    figure(10) % mean
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(dist_list/1000000,dist_mean_BCWLS/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
    xlabel('Distance (km $\times$ 10$^3$) ','Interpreter','latex')
    ylabel('Mean absolute range error (km)')
    box on
    
    figure(11) % standard deviation
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(dist_list/1000000,dist_std_BCWLS/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
    xlabel('Distance (km $\times$ 10$^3$)','Interpreter','latex')
    ylabel('STD absolute range error (km)')
    box on
    
    figure(12) % Coordinate position errors
    subplot(3,1,1); hold on; plot(dist_list/1000000,pos_mean_LOS/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
    subplot(3,1,2); hold on; plot(dist_list/1000000,pos_mean_HCA/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
    subplot(3,1,3); hold on; plot(dist_list/1000000,pos_mean_VCA/1000,'Color',shade,'LineStyle','-.','LineWidth',1.5)
    
    solar_idx = solar_idx+1; % update index value
    end

    figure(12)
    subplot(3,1,1); box on; ax = gca; ax.FontSize = 14; ax.FontName = 'Times New Roman';
    subplot(3,1,2); box on; ax = gca; ax.FontSize = 14; ax.FontName = 'Times New Roman';
    subplot(3,1,3); box on; ax = gca; ax.FontSize = 14; ax.FontName = 'Times New Roman';
    han=axes(figure(12),'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Absolute position error (km)');
    xlabel(han,'Distance (km $\times$ 10$^3$)','Interpreter','latex');
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
end

%% Sorting Experiment

if test_sort

    fprintf("------------------------- Sorting Test -------------------------\n")
    
    distance_sc = distance_sorting*10^3;           % (m) initialize
    r_sc_I = inertial_uv*(distance_sc); % (m) Inertial position vector

    % Get detected crater uncertainty
    [~,repeat_matrix_detections,~,~,~] = angular_error_calc(distance_sc/1000, moon_angle);
    mat_size = size(repeat_matrix_detections);
    num_craters = mat_size(1); 

    % Initialize errors 
    dist_LS_errors_std = zeros(num_craters, N); % angle std

    dist_WLS_errors_std = zeros(num_craters, N);

    dist_BCWLS_errors_std = zeros(num_craters, N);

    dist_LS_errors_rad = zeros(num_craters, N); % crater radius

    dist_WLS_errors_rad = zeros(num_craters, N);

    dist_BCWLS_errors_rad = zeros(num_craters, N);

    dist_LS_errors_az = zeros(num_craters, N); % measured azimuth

    dist_WLS_errors_az = zeros(num_craters, N);

    dist_BCWLS_errors_az = zeros(num_craters, N);

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

    fprintf("%.f detected craters\n\n", num_craters)

    for idx = 3:num_craters % use all craters

        fprintf("Computing %.f best craters\n", idx)

        for i = 1:N % Monte Carlo iterations

        % Generate angular measurements and get uncertainty
        [std_elevation_og,std_azimuth_og,meas_azimuth_og,meas_elevation_og,r_crater_LOS_nom_og,R_LOS_I] = crater_pos(r_sc_I,sc_bearings,std_dev_unnorm);

        % Sort values by angle std
        [~, idx_sort_std] = sort(std_dev_unnorm,'ascend');
        meas_azimuth = meas_azimuth_og(idx_sort_std);
        meas_elevation = meas_elevation_og(idx_sort_std);
        r_crater_LOS_nom = r_crater_LOS_nom_og(:,idx_sort_std);                                  
        std_azimuth = std_azimuth_og(idx_sort_std);
        std_elevation = std_elevation_og(idx_sort_std);

        meas_azimuth = meas_azimuth(1:idx); % get idx best values
        meas_elevation = meas_elevation(1:idx);
        r_crater_LOS_nom = r_crater_LOS_nom(:,1:idx);
        std_azimuth = std_azimuth(1:idx);
        std_elevation = std_elevation(1:idx);
        

        % LS estimate error
        [pos_LS_est,LS_error,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation,r_crater_LOS_nom); 

        % WLS estimate error
        [pos_WLS_est,WLS_error,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,r_crater_LOS_nom,pos_LS_est,A_mat,z_mat);

        % BCWLS esimtate error
        [~,BCWLS_error] = BCWLS_estimate(std_azimuth,std_elevation,meas_azimuth,meas_elevation,pos_WLS_est,A_mat,W_block,do);

        % Store values
        dist_LS_errors_std(idx,i) = LS_error;       % (m) Range errors
        dist_WLS_errors_std(idx,i) = WLS_error;
        dist_BCWLS_errors_std(idx,i) = BCWLS_error;


        % Sort values by crater radius
        [~, idx_sort_rad] = sort(crater_radius,'descend');
        meas_azimuth = meas_azimuth_og(idx_sort_rad);
        meas_elevation = meas_elevation_og(idx_sort_rad);
        r_crater_LOS_nom = r_crater_LOS_nom_og(:,idx_sort_rad);                                  
        std_azimuth = std_azimuth_og(idx_sort_rad);
        std_elevation = std_elevation_og(idx_sort_rad);

        meas_azimuth = meas_azimuth(1:idx); % get idx best values
        meas_elevation = meas_elevation(1:idx);
        r_crater_LOS_nom = r_crater_LOS_nom(:,1:idx);
        std_azimuth = std_azimuth(1:idx);
        std_elevation = std_elevation(1:idx);
        

        % LS estimate error
        [pos_LS_est,LS_error,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation,r_crater_LOS_nom); 

        % WLS estimate error
        [pos_WLS_est,WLS_error,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,r_crater_LOS_nom,pos_LS_est,A_mat,z_mat);

        % BCWLS esimtate error
        [~,BCWLS_error] = BCWLS_estimate(std_azimuth,std_elevation,meas_azimuth,meas_elevation,pos_WLS_est,A_mat,W_block,do);

        % Store values
        dist_LS_errors_rad(idx,i) = LS_error;       % (m) Range errors
        dist_WLS_errors_rad(idx,i) = WLS_error;
        dist_BCWLS_errors_rad(idx,i) = BCWLS_error;


        % Sort values by measured azimuth
        [~, idx_sort_az] = sort(meas_azimuth,'ascend');
        meas_azimuth = meas_azimuth_og(idx_sort_az);
        meas_elevation = meas_elevation_og(idx_sort_az);
        r_crater_LOS_nom = r_crater_LOS_nom_og(:,idx_sort_az);                                  
        std_azimuth = std_azimuth_og(idx_sort_az);
        std_elevation = std_elevation_og(idx_sort_az);

        meas_azimuth = meas_azimuth(1:idx); % get idx best values
        meas_elevation = meas_elevation(1:idx);
        r_crater_LOS_nom = r_crater_LOS_nom(:,1:idx);
        std_azimuth = std_azimuth(1:idx);
        std_elevation = std_elevation(1:idx);
        

        % LS estimate error
        [pos_LS_est,LS_error,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation,r_crater_LOS_nom); 

        % WLS estimate error
        [pos_WLS_est,WLS_error,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,r_crater_LOS_nom,pos_LS_est,A_mat,z_mat);

        % BCWLS esimtate error
        [~,BCWLS_error] = BCWLS_estimate(std_azimuth,std_elevation,meas_azimuth,meas_elevation,pos_WLS_est,A_mat,W_block,do);

        % Store values
        dist_LS_errors_az(idx,i) = LS_error;       % (m) Range errors
        dist_WLS_errors_az(idx,i) = WLS_error;
        dist_BCWLS_errors_az(idx,i) = BCWLS_error;
        end

    end

    % Calculate Mean values
    dist_mean_LS_std = mean(dist_LS_errors_std,2);         % LS
    dist_mean_LS_rad = mean(dist_LS_errors_rad,2);
    dist_mean_LS_az = mean(dist_LS_errors_az,2);

    dist_mean_WLS_std = mean(dist_WLS_errors_std,2);       % WLS
    dist_mean_WLS_rad = mean(dist_WLS_errors_rad,2);
    dist_mean_WLS_az = mean(dist_WLS_errors_az,2);

    dist_mean_BCWLS_std = mean(dist_BCWLS_errors_std,2);   % BCWLS
    dist_mean_BCWLS_rad = mean(dist_BCWLS_errors_rad,2);
    dist_mean_BCWLS_az = mean(dist_BCWLS_errors_az,2);

    % Plots

    figure % std
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(3:num_craters,dist_mean_LS_std(3:end)/1000,'r-','LineWidth', 1.5)
    plot(3:num_craters,dist_mean_WLS_std(3:end)/1000,'k-','LineWidth', 1.5)
    plot(3:num_craters,dist_mean_BCWLS_std(3:end)/1000,'g--','LineWidth', 1.5)
    xlabel('Number of craters')
    ylabel('Mean absolute range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',14)
    box on
    
    figure % crater radius
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(3:num_craters,dist_mean_LS_rad(3:end)/1000,'r-','LineWidth', 1.5)
    plot(3:num_craters,dist_mean_WLS_rad(3:end)/1000,'k-','LineWidth', 1.5)
    plot(3:num_craters,dist_mean_BCWLS_rad(3:end)/1000,'g--','LineWidth', 1.5)
    xlabel('Number of craters')
    ylabel('Mean absoulte range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',14)
    box on

    figure % azimuth
    hold on
    ax = gca;
    ax.FontSize = 14;
    ax.FontName = 'Times New Roman';
    plot(3:num_craters,dist_mean_LS_az(3:end)/1000,'r-','LineWidth', 1.5)
    plot(3:num_craters,dist_mean_WLS_az(3:end)/1000,'k-','LineWidth', 1.5)
    plot(3:num_craters,dist_mean_BCWLS_az(3:end)/1000,'g--','LineWidth', 1.5)
    xlabel('Number of craters')
    ylabel('Mean absolute range error (km)')
    legend('LS','WLS','BCWLS', 'FontName','Times New Roman', 'FontSize',14)
    box on
end

