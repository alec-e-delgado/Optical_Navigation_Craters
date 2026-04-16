function [std_elevation,std_azimuth,azimuth_meas,elevation_meas, r_crater_LOS_nom,R_auxLOS_I] = ...
    crater_pos(r_sc_I, sc_bearing, std_expanded)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRATER_POS
%  Calculaes the positions of detected craters defined in an inertial 
%  reference frame centered about the Moon as a function of the altitude of 
%  the space craft and the incidence angle between the detcted crater and
%  spacecraft, azimuth and elevation measurements with simulated noise are
%  then generated
%
% INPUTS
%   r_sc_I           - (m) inertial position of the spacecraft
%   sc_bearing       - (deg) bearing angles of spacecraft and craters 
%   std_expanded     - (rad) standard deviation of angular error for each crater 
%
% OUTPUTS
%   std_elevation    - (rad) analytical standard deviation of elevation angle
%   std_azimuth      - (rad) analytical standard deviation of azimuth angle
%   azimuth_meas     - (rad) measured azimuth angle in LOS frame
%   elevation_meas   - (rad) measured elevation angle in LOS frame
%   r_crater_LOS_nom - (m)   known location of detected craters
%   R_aux_LOS_I      - rotation matrix from LOS frame to Inertial frame
%
% ADDITIONAL INFORMATION
%  - theta: a randomly generated angle from the (+) y_aux axis that locates
%  the crater along a given circle of possibilites
%  - auxiliary refernce frame: An intermittent frame used to calculate the
%  crater location
%    * x_aux: points from center of moon to sc position
%    * y_aux: some direction normal to z_I and x_aux
%    * z_aux: some direction normal tp x_aux and z_aux
%  - inertial reference frame: a non-rotating reference frame centered about
%  the moon
%  - Note that this script locates craters using a known location for the
%  spacecraft (unit vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find nominal alpha that gives true position of crater

% Define Parameters
radius_Moon = 1.7374e6;             % (m)
sc_bearing_rad = sc_bearing*pi/180; % (rad)

% Calculate angle between r_sc_I and r_sc_I - r_c_I (alpha) using law of sines
distance_sc = norm(r_sc_I);                                       % (m) distance from center of moon
alpha_nom = asin(radius_Moon*sin(pi-sc_bearing_rad)/distance_sc); % (rad) nominal

% Calculate angle between r_sc_I and r_c_I (gamma) 
beta_nom = pi - sc_bearing_rad;        % (rad) angle between r_c_I and r_sc_I - r_c_I
gamma_nom = pi - beta_nom - alpha_nom; % (rad)

%% Calculate crater position in Auxiliary reference frame
% Aux reference frame has its origin at the center of the moon with its x-axis 
% pointing towards the spacecraft

% Random angles to place crater positions along a circle
n = length(sc_bearing_rad);    % number of elements
theta_crater = 2*pi*rand(n,1); % (rad) angles for crater position generation

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

% Calculate angular errors 
b_ray = std_expanded/sqrt((4-pi)/2); 
angular_errors = raylrnd(b_ray);                   % (rad) rayleigh distribution
rand_theta = 2*pi* rand(length(angular_errors),1); % rand angles for error generation

%% Calculate perturbed crater positions in err reference frame 
% ERR reference frame has its orgin centered at the spacecraft with its x
% axis pointing towards the nominal crater position

% Pre allocate
r_Ecrater_LOS = zeros(3, length(r_crater_LOS_nom)); % initialize

for i= 1:length(r_crater_LOS_nom)
    r_sc_c_dist = norm(r_crater_LOS_nom(:,i)); % (m) from spacecraft to crater
    
    % Calculate position of perturbed craters (Ecrater) in ERR frame
    x_Ecrater_ERR = cos(angular_errors(i));                             
    y_Ecrater_ERR = sin(angular_errors(i)) .* cos(rand_theta(i));
    z_Ecrater_ERR = sin(angular_errors(i)) .* sin(rand_theta(i));
    
    r_Ecrater_ERR = r_sc_c_dist*[x_Ecrater_ERR y_Ecrater_ERR z_Ecrater_ERR]; % extend by nominal distance
    
    % Generate rotation matrix using unit vectors defined in LOS frame
    x_uv_ERR = r_crater_LOS_nom(:,i)/norm(r_crater_LOS_nom(:,i));
    y_uv_ERR = cross([0 0 1], x_uv_ERR)/norm(cross([0 0 1], x_uv_ERR));
    z_uv_ERR = cross(x_uv_ERR, y_uv_ERR);
    
    R_ERR_LOS = [x_uv_ERR, y_uv_ERR', z_uv_ERR']; % rotate from ERR to LOS frames
    
    % Ecrater positions in LOS frame
    r_Ecrater_LOS(:,i) = R_ERR_LOS*r_Ecrater_ERR';
end

% Ecrater positions in I frame
r_Ecrater_I = R_auxLOS_I*(r_Ecrater_LOS-[0;distance_sc;0]);

%% Calculate Deviations in Azimuth and Elevation angles
% These angles are found in the LOS reference frame

% Calculate azimuth and elevation angles of nominal crater position
azimuth_nom = atan2(r_crater_LOS_nom(2,:), r_crater_LOS_nom(1,:)); % (rad) nominal angle
elevation_nom = atan2(r_crater_LOS_nom(3,:), sqrt(r_crater_LOS_nom(1,:).^2+r_crater_LOS_nom(2,:).^2)); % (rad) nominal angle

% Calculate azimuth and elevation angles of Ecraters
azimuth_meas = atan2(r_Ecrater_LOS(2,:), r_Ecrater_LOS(1,:)); % (rad) measured angle
elevation_meas = atan2(r_Ecrater_LOS(3,:), sqrt(r_Ecrater_LOS(1,:).^2+r_Ecrater_LOS(2,:).^2)); % (m) measured angle

% Analytical Standard Deviations
std_elevation = std_expanded/sqrt((4-pi)/2); % (rad) from paper
std_azimuth = std_expanded./sqrt((4-pi)/2)./cos(elevation_nom'); % (rad) from paper

end
