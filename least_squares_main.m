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
for j = 1:rep_matrix_size(1)
    fprintf('%d -- sc_inc = %.2f -- crater_radius = %.2f\n', j, repeat_matrix_detections(j,3),repeat_matrix_detections(j,4))
end

desired_indices = 1:46; % list of indices to use for crater generation

% Extract crater quantities associated with desired indeces
std_devs = repeat_matrix_detections(desired_indices,2); % rad
means = repeat_matrix_detections(desired_indices,1); % rad
sc_bearings = repeat_matrix_detections(desired_indices, 3); % rad

% Pass statistic values 
[r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearings, std_devs, means);