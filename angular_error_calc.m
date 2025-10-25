function [angular_errors, x_2, y_2, z_2] = angular_error_calc(distance, moon_angle)
%ANGULAR_ERROR_CALC generates angular errors along a normal distribution
%
% INPUTS
%   - distance: the altitude above the lunar surface (m)
%   - moon_angle: the incidence angle with the lunar surface (degrees)
%
% OUTPUTS
%   - angular_errors: angular errors generated from the statistics
%   associated with the detected craters following a normal distribution
%   -x_2,y_2,z_2: coordinates of vectors that display the angular errors
%   from the true angles of the craters
%
% ADDITIONAL INFORMATION
%   - This function loads the file 'discretized_results.mat' that provides
%   the data and statistics of the detected craters at various altitudes
%   and moon angles 
%

    load('discretized_results.mat')
    
    % --- Select indices where range ≤ distance and moon-angle ≤ moon_angle ---
    vec_find_range        = find(discretize_range <= distance);
    vec_find_moonangle    = find(discretize_moon_angle <= moon_angle);
    
    % --- Slice the 8-D dataset at those index limits ---
    sliced_data = discretized_data(2,2,1:3, ...
        vec_find_moonangle(end), vec_find_range(end), :,:,2);
    
    % --- Extract the detection counts matrix ---
    num_data = squeeze(sliced_data(1,1,1,1,1,:,:,1)); 
    
    % --- Find non-zero entries and their linear indices ---
    [vec_find_nonzero_row, vec_find_nonzero_col] = find(num_data ~= 0);
    linear_idx = sub2ind(size(num_data), vec_find_nonzero_row, vec_find_nonzero_col);
    
    % --- Pull out mean/std matrices and select only non-zero points ---
    mean_data_matrix = squeeze(sliced_data(1,1,2,1,1,:,: ,1));
    std_data_matrix  = squeeze(sliced_data(1,1,3,1,1,:,: ,1));
    mean_data        = mean_data_matrix(linear_idx);
    std_data         = std_data_matrix(linear_idx);
    number_nonzero   = num_data(linear_idx); % each element is the number of detected craters with specific radius and incidence angle
    
    % --- Convert matrix indices into physical angles & radii ---
    sc_inc_nonzero = discretize_sc_inc_angle(vec_find_nonzero_row);
    radius_nonzero = discretize_radius(vec_find_nonzero_col); % elements are mapped to number_nonzero
    
    % --- Round detection counts to integer crater numbers ---
    num_craters = round(number_nonzero);
    
    % --- Repeat each statistic by crater for sampling ---
    expanded_mean_data = repelem(mean_data(:), num_craters); % each element now corresponds to one crater
    expanded_std_data  = repelem(std_data(:),  num_craters);
    expanded_sc_inc    = repelem(sc_inc_nonzero(:), num_craters);
    expanded_radius    = repelem(radius_nonzero(:), num_craters);
    repeat_matrix_detections = [expanded_mean_data, expanded_std_data, expanded_sc_inc, expanded_radius];
    
    % --- Generate angular error deviations for each crater ---
    angular_errors = normrnd(repeat_matrix_detections(:,1), repeat_matrix_detections(:,2));
    
    % --- Create random angles [0, 2π) ---
    rand_theta = 2 * pi() * rand(1, length(angular_errors));
    rand_theta = rand_theta';
    angles_matrix = [angular_errors, rand_theta];
    
    % --- Project spherical offsets into 3D vectors ---
    x_2 = cos(angles_matrix(:,1));                             
    y_2 = sin(angles_matrix(:,1)) .* cos(angles_matrix(:,2));
    z_2 = sin(angles_matrix(:,1)) .* sin(angles_matrix(:,2));
end 
