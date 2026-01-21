function [r_crater_I_nom, r_crater_aux_nom] = crater_pos(r_sc_I, sc_bearing, std_expanded, mean_expanded)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRATER_POS
%  - This script calculaes the positions of detected craters defined in
%  an inertial reference frame centered about the Moon as a function of the
%  altitude of the space craft and the incidence angle between the detcted
%  crater and spacecraft
%
% INPUTS
%  - altitude: the altitude above the lunar surface (m)
%  - sc_bearing: bearing angles of spacecraft and crater (deg) 
%  - std_expanded: standard deviation of angular error for each crater (rad)
%  - mean_expanded: mean of engular error for each crater (rad)
%
% OUTPUTS
%  - r_crater_I: location of crater in intertial reference frame (m)
%
% ADDITIONAL INFORMATION
% - theta: a randomly generated angle from the (+) y_aux axis that locates
% the crater along a given circle of possibilites
% - auxiliary refernce frame: An intermittent frame used to calculate the
% crater location
%   * x_aux: points from center of moon to sc position
%   * y_aux: some direction normal to z_I and x_aux
%   * z_aux: some direction normal tp x_aux and z_aux
% - inertial reference frame: a non-rotating reference frame centered about
% the moon
% - Note that this script locates craters using a known location for the
% spacecraft (unit vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define seed value for consistent theta
% rng(3);

%% Find nominal alpha that gives true position of crater

% Define Parameters
radius_Moon = 1.7374e6; % m
sc_bearing_rad = sc_bearing*pi/180; % convert to rad

% Calculate angle between r_sc_I and r_sc_I - r_c_I (alpha) using law of sines
distance_sc = norm(r_sc_I); % distance from center of moon
alpha_nom = asin(radius_Moon*sin(pi-sc_bearing_rad)/distance_sc);

% Calculate angle between r_sc_I and r_c_I (gamma) 
beta_nom = pi - sc_bearing_rad; % angle between r_c_I and r_sc_I - r_c_I
gamma_nom = pi - beta_nom - alpha_nom; % (rad)

%% Calculate crater position in Auxiliary reference frame
% Aux reference frame has its origin at the center of the moon with its x-axis 
% pointing towards the spacecraft

% Random angles to place crater positions along a circle
n = length(sc_bearing_rad); % number of elements
theta_crater = 2*pi*rand(n,1); % (rad) get random angles

% Crater positoin in Aux
x_crater_aux_nom = cos(gamma_nom)*radius_Moon; % distance along vector pointing to sc
y_crater_aux_nom = sin(gamma_nom).*cos(theta_crater)*radius_Moon; 
z_crater_aux_nom = sin(gamma_nom).*sin(theta_crater)*radius_Moon;

r_crater_aux_nom = [x_crater_aux_nom y_crater_aux_nom z_crater_aux_nom]';

% Generate rotation matrix using unit vectors defined in inertial refernece frame 
x_uv_aux = r_sc_I/distance_sc;
y_uv_aux = cross([0 0 1], x_uv_aux)/norm(cross([0 0 1], x_uv_aux)); 
z_uv_aux = cross(x_uv_aux, y_uv_aux);

R_aux_I = [x_uv_aux, y_uv_aux', z_uv_aux']; % Rotation matrix from aux to I

% Crater position in inertial
r_crater_I_nom = R_aux_I*r_crater_aux_nom; 

%% Calculate crater position in LOS reference frame
% LOS reference frame has its origin centered at the spacecraft with its y
% axis pointing towards center of the moon

% Crater position in LOS
y_crater_LOS_nom = distance_sc-x_crater_aux_nom;
x_crater_LOS_nom = y_crater_aux_nom;
z_crater_LOS_nom = z_crater_aux_nom;

r_crater_LOS_nom = [x_crater_LOS_nom y_crater_LOS_nom z_crater_LOS_nom]';

% Generate rotation matrix using unit vectors defined in inertial reference frame
y_uv_LOS = -x_uv_aux;
x_uv_LOS = y_uv_aux;
z_uv_LOS = z_uv_aux;

R_auxLOS_I = [x_uv_LOS', y_uv_LOS, z_uv_LOS']; % Rotation matrix from LOS to I

% Crater position in Inretial used for validation
r_crater_I2_nom = R_auxLOS_I*(r_crater_LOS_nom-[0;distance_sc;0]);
% 
% % Used to visualize the aux and LOS reference frames
% x_origin = 0;
% y_origin = 0;
% z_origin = 0;
% figure; % Open a new figure window
% hold on 
% for i = 1:3
%     if i == 1
%         u = x_uv_aux(1);
%         v = x_uv_aux(2);
%         w = x_uv_aux(3);
%         % Plot the vector
%         quiver3(x_origin, y_origin, z_origin, u, v, w, 0,'Color', 'r'); % The '0' disables automatic scaling
%     elseif i == 2
%         u = y_uv_aux(1);
%         v = y_uv_aux(2);
%         w = y_uv_aux(3);
%         % Plot the vector
%         quiver3(x_origin, y_origin, z_origin, u, v, w, 0,'Color', 'g'); % The '0' disables automatic scaling
%     else
%         u = z_uv_aux(1);
%         v = z_uv_aux(2);
%         w = z_uv_aux(3);
%         % Plot the vector
%         quiver3(x_origin, y_origin, z_origin, u, v, w, 0,'Color', 'b'); % The '0' disables automatic scaling
%     end
% end
% axis equal; % Use equal data unit lengths along each axis
% xlim([-1 1]); % Set x-axis limits
% ylim([-1 1]); % Set y-axis limits
% zlim([-1 1]); % Set z-axis limits
% title('3D Unit Vector Plot');
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% grid on;
% 
% figure 
% hold on
% for i = 1:3
%     if i == 1
%         u = x_uv_LOS(1);
%         v = x_uv_LOS(2);
%         w = x_uv_LOS(3);
%         % Plot the vector
%         quiver3(x_origin, y_origin, z_origin, u, v, w, 0,'Color', 'm'); % The '0' disables automatic scaling
%     elseif i == 2
%         u = y_uv_LOS(1);
%         v = y_uv_LOS(2);
%         w = y_uv_LOS(3);
%         % Plot the vector
%         quiver3(x_origin, y_origin, z_origin, u, v, w, 0,'Color', 'y'); % The '0' disables automatic scaling
%     else
%         u = z_uv_LOS(1);
%         v = z_uv_LOS(2);
%         w = z_uv_LOS(3);
%         % Plot the vector
%         quiver3(x_origin, y_origin, z_origin, u, v, w, 0,'Color', 'c'); % The '0' disables automatic scaling
%     end
% end
% axis equal; % Use equal data unit lengths along each axis
% xlim([-1 1]); % Set x-axis limits
% ylim([-1 1]); % Set y-axis limits
% zlim([-1 1]); % Set z-axis limits
% title('3D Unit Vector Plot');
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% grid on;


% Calculate angular errors 
% angular_errors = normrnd(0, std_expanded);
% angular_errors = normrnd(mean_expanded, std_expanded);
b_ray = std_expanded/sqrt((4-pi)/2);
angular_errors = raylrnd(b_ray);
rand_theta = 2*pi* rand(1,length(angular_errors)); % random angles in circle
rand_theta = rand_theta';

%% Calculate perturbed crater positions in err reference frame 
% ERR reference frame has its orgin centered at the spacecraft with its x
% axis pointing towards the nominal crater position

r_sc_c_dist = norm(r_crater_LOS_nom); % distance from spacecraft to crater

% Calculate position of perturbed craters (Ecrater) in ERR frame
x_Ecrater_ERR = cos(angular_errors);                             
y_Ecrater_ERR = sin(angular_errors) .* cos(rand_theta);
z_Ecrater_ERR = sin(angular_errors) .* sin(rand_theta);

r_Ecrater_ERR = r_sc_c_dist*[x_Ecrater_ERR y_Ecrater_ERR z_Ecrater_ERR]; % extend ny nominal distance

% Generate rotation matrix using unit vectors defined in LOS frame
x_uv_ERR = r_crater_LOS_nom/norm(r_crater_LOS_nom);
y_uv_ERR = cross([0 0 1], x_uv_ERR)/norm(cross([0 0 1], x_uv_ERR));
z_uv_ERR = cross(x_uv_ERR, y_uv_ERR);

R_ERR_LOS = [x_uv_ERR, y_uv_ERR', z_uv_ERR']; % rotate from ERR to LOS frames

% Ecrater positions in LOS frame
r_Ecrater_LOS = R_ERR_LOS*r_Ecrater_ERR';

% Ecrater positions in I frame
r_Ecrater_I = R_auxLOS_I*(r_Ecrater_LOS-[0;distance_sc;0]);

%% Calculate Deviations in Azimuth and Elevation angles
% These angles are found in the LOS reference frame

% Calculate azimuth and elevation angles of nominal crater position
azimuth_nom = atan2(r_crater_LOS_nom(2), r_crater_LOS_nom(1));
elevation_nom = atan2(r_crater_LOS_nom(3), sqrt(r_crater_LOS_nom(1)^2+r_crater_LOS_nom(2)^2));

% Calculate azimuth and elevation angles of Ecraters
azimuth = atan2(r_Ecrater_LOS(2,:), r_Ecrater_LOS(1,:));
elevation = atan2(r_Ecrater_LOS(3,:), sqrt(r_Ecrater_LOS(1,:).^2+r_Ecrater_LOS(2,:).^2));

mean_azimuth_num = mean(azimuth);
mean_elevation_num = mean(elevation);

% Numerical Standard Deviations 
std_azimuth_num = std(azimuth);
std_elevation_num = std(elevation);

% Analytical Standard Deviations
% std_azimuth_analytical = (std_expanded(1)/sqrt(2))/cos(azimuth_nom);
% std_elevation_analytical = std_expanded(1)/sqrt(2)/cos(elevation_nom);
std_elevation_analytical = std_expanded(1)/sqrt((4-pi)/2);
std_azimuth_analytical = std_expanded(1)/sqrt((4-pi)/2)/cos(elevation_nom);

num_error = length(std_expanded);

% figure % used to visualize crater error locations in error crater frame
% % quiver3(zeros(num_error,1), zeros(num_error,1), zeros(num_error,1), x_error_aux_uv, y_error_aux_uv, z_error_aux_uv, 'r');
% quiver3(zeros(num_error,1), zeros(num_error,1), zeros(num_error,1), r_Ecrater_ERR(:,1), ...
%      r_Ecrater_ERR(:,2),r_Ecrater_ERR(:,3), 'r');
% xlabel('x-axis')
% ylabel('y-axis')
% zlabel('z-axis')
% axis equal


theta_circle = 0:0.01:2*pi; 
x_circle_aux = cos(gamma_nom(1))*radius_Moon*ones(length(theta_circle), 1); % used to plot circle of possibilites
y_circle_aux = sin(gamma_nom(1))*cos(theta_circle)*radius_Moon; 
z_circle_aux = sin(gamma_nom(1))*sin(theta_circle)*radius_Moon;
r_circle_aux = [x_circle_aux'; y_circle_aux; z_circle_aux];
r_circle_I = R_aux_I*r_circle_aux;

[xsphere, ysphere, zsphere] = sphere(30);
xsphere = radius_Moon * xsphere;
ysphere = radius_Moon * ysphere;
zsphere = radius_Moon * zsphere;

% figure;
% axis equal
% moon_sphere = surf(xsphere, ysphere, zsphere);
% set(moon_sphere, 'FaceColor', [0.5 0.5 0.5]); 
% hold on 
% plot3(x_crater_aux2_nom, y_crater_aux2_nom, z_crater_aux2_nom, 'r*', 'MarkerSize',5, 'LineWidth',1.5)
% plot3(x_circle_aux, y_circle_aux, z_circle_aux, 'm-', 'MarkerSize',5, 'LineWidth',2)
% plot3(distance_sc, 0, 0, 'b*', 'MarkerSize',5, 'LineWidth',2)
% xlabel('x-axis (m)')
% ylabel('y-axis (m)')
% zlabel('z-axis (m)')
% title('Auxiliary reference frame')
% text(distance_sc, 0, 0, 'SC', 'Color',[1 1 1])

% % plot in the inertial reference frame
% figure;
% axis equal
% moon_sphere = surf(xsphere, ysphere, zsphere);
% set(moon_sphere, 'FaceColor', [0.5 0.5 0.5]); 
% hold on 
% % plot3(r_crater_I_nom(1,:), r_crater_I_nom(2,:), r_crater_I_nom(3,:), 'yx', 'MarkerSize',5, 'LineWidth',1.5)
% plot3(r_Ecrater_I(1,:), r_Ecrater_I(2,:), r_Ecrater_I(3,:), 'yx', 'MarkerSize',5, 'LineWidth',1.5)
% plot3(r_circle_I(1,:), r_circle_I(2,:), r_circle_I(3,:), 'g-', 'MarkerSize',5, 'LineWidth',1.5)
% plot3(r_sc_I(1),r_sc_I(2), r_sc_I(3), 'b*', 'MarkerSize',5, 'LineWidth',2)
% xlabel('x-axis (m)')
% ylabel('y-axis (m)')
% zlabel('z-axis (m)')
% title('Inertial reference frame centered at the Moon')
% text(r_sc_I(1), r_sc_I(2), r_sc_I(3), 'SC', 'Color',[1 1 1])

fprintf("Nominal crater position (LOS) (%.2f, %.2f, %.2f) km \n", r_crater_LOS_nom(1)/1000, r_crater_LOS_nom(2)/1000, r_crater_LOS_nom(3)/1000)
fprintf("Bearing angle: %.2f deg \n \n", sc_bearing)
fprintf("Nominal Azimuth: %.4f deg ---- Nominal Elevation: %.4f deg \n", rad2deg(azimuth_nom), rad2deg(elevation_nom))
fprintf("Mean Azimuth:    %.4f deg ---- Mean Elevatiopn:   %.4f deg \n \n", rad2deg(mean_azimuth_num), rad2deg(mean_elevation_num))
fprintf("Numerical Azimuth std:  %.4f deg ---- Numerical Elevation std:  %.4f deg \n", rad2deg(std_azimuth_num), rad2deg(std_elevation_num))
fprintf("Analytical Azimuth std: %.4f deg ---- Analytical Elevation std: %.4f deg \n\n", rad2deg(std_azimuth_analytical), rad2deg(std_elevation_analytical))
fprintf("Angular error Mean:  %.4f deg ---- Angular Error std:  %.4f deg \n", rad2deg(mean_expanded(1)), rad2deg(std_expanded(1)))


end
