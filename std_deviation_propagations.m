close all; clear; clc;

% For numerically finding the sigmas, we will fix the sc_bearing of the
% crater, the crater radius (20 m), and theta in crater_pos such that the
% same crater is being considered. With these fixed values, the values for
% the mean and standard deviation of the angular error are acquired from
% angular_error_calc (just 

% Define parameters
altitude_sc = 42000; % current distance away from moon surface (need units)
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

% Used to select desired index

% % Display list of craters
% for j = 1:rep_matrix_size(1)
%     fprintf('%d -- sc_inc = %.2f -- crater_radius = %.2f\n', j, repeat_matrix_detections(j,3),repeat_matrix_detections(j,4))
% end

% % Select which crater to propagate numerically
% while is_pass == false
%     desired_index = input('Enter desired index: ');
%     if desired_index > 0 && desired_index <= rep_matrix_size(1)
%         is_pass = true;
%     else
%         fprintf('Invalid index. Max size is %d Please try again.\n', rep_matrix_size(1));
%     end
% end

des_idx = 1:num_craters;

% Extract crater quantities associated with desired index
std_dev = repeat_matrix_detections(des_idx,2); % rad
mean = repeat_matrix_detections(des_idx,1); % rad
sc_bearing = repeat_matrix_detections(des_idx, 3); % rad

N = 10000; % number of runs in monte carlos simulation

num_sim = 30; % Number of monte carlos simulations


% Info to store across runs
stack_mat = [];


for j = 1:num_sim
    temp_mat = zeros(num_craters,15);

    for i = 1:num_craters
    
        std_dev_expanded = ones(N,1)*std_dev(i);
        mean_expanded = ones(N,1)*mean(i);
    
        [info_store, r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearing(i), std_dev_expanded, mean_expanded);

        temp_mat(i,1:12) = [i, info_store];
    end
    stack_mat = cat(1,stack_mat,temp_mat);
    stack_mat(:,13) = abs(stack_mat(:,7)-stack_mat(:,9))./stack_mat(:,7)*100;
    stack_mat(:,14) = abs(stack_mat(:,8)-stack_mat(:,10))./stack_mat(:,8)*100;
    stack_mat(:,15) = stack_mat(:,12)./stack_mat(:,11);
    % stack_mat(end+1,:) = zeros(1,15);
end

% Make matrices for each crater

heading = {'Crater', 'Bearing Angle',' Aximuth nom', 'Elevation Nom', 'Mean Azimuth', 'Mean Elevation',...
           'Azimuth std num', 'Elevation std num','Azimuth std Analytical','Elevation std Analytical',...
           'Angle Mean','Angle std','Azimuth std % error','Elevation std % error','Ratio'};

filename = 'testdata.xlsx';
writecell(heading,filename,'Sheet',1)
writematrix(stack_mat,filename,'Sheet',1,'Range','A2')

% Organize by craters
crater_mat = sortrows(stack_mat,1,'ascend');
writecell(heading,filename,'Sheet',2)
writematrix(crater_mat,filename,'Sheet',2,'Range','A2')

% Organize by bearing angle
angle_mat = sortrows(stack_mat,2,'ascend');
writecell(heading,filename,'Sheet',3)
writematrix(angle_mat,filename,'Sheet',3,'Range','A2')

% Organize by elevation
elevation_mat = sortrows(stack_mat,4,'descend','ComparisonMethod','abs');
writecell(heading,filename,'Sheet',4)
writematrix(elevation_mat,filename,'Sheet',4,'Range','A2')





