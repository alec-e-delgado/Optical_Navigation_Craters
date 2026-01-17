close all; clear; clc;

% For numerically finding the sigmas, we will fix the sc_bearing of the
% crater, the crater radius (20 m), and theta in crater_pos such that the
% same crater is being considered. With these fixed values, the values for
% the mean and standard deviation of the angular error are acquired from
% angular_error_calc (just 

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

% Select which crater to propagate numerically
while is_pass == false
    desired_index = input('Enter desired index: ');
    if desired_index > 0 && desired_index <= rep_matrix_size(1)
        is_pass = true;
    else
        fprintf('Invalid index. Max size is %d Please try again.\n', rep_matrix_size(1));
    end
end

% Extract crater quantities associated with desired index
std_dev = repeat_matrix_detections(desired_index,2); % rad
mean = repeat_matrix_detections(desired_index,1); % rad
sc_bearing = repeat_matrix_detections(desired_index, 3); % rad

% Expand parameters for number of desired samples
N = 1000000; % number of samples to be taken
std_dev_expanded = ones(N,1)*std_dev;
mean_expanded = ones(N,1)*mean;

[r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearing, std_dev_expanded, mean_expanded);