function [r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearing)
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

% Define Parameters
radius_Moon = 1.7374e6; % m
sc_bearing_rad = sc_bearing*pi/180; % convert to rad

% Calculate angle between r_sc_I and r_sc_I - r_c_I (alpha) using law of sines
distance_sc = norm(r_sc_I); % magnitude
alpha = asin(radius_Moon*sin(pi-sc_bearing_rad)/distance_sc);

% Calculate angle between r_sc_I and r_c_I (gamma) 
beta = 180*pi/180 - sc_bearing_rad; % angle between r_c_I and r_sc_I - r_c_I
gamma = pi - beta - alpha; % (rad)

% Define auxilary reference frame where gamma is the angle between the sc
% and crater from the center of the Moon
n = length(sc_bearing_rad); % number of elements
theta_rand = 2*pi*rand(n,1); % (rad) get random angles

x_crater_aux = cos(gamma)*radius_Moon; % distance along vector pointing to sc
y_crater_aux = sin(gamma).*cos(theta_rand)*radius_Moon; 
z_crater_aux = sin(gamma).*sin(theta_rand)*radius_Moon;

% Generate rotation matrix using unit vectors for auxiliary axes 
x_uv_aux = r_sc_I/distance_sc;
y_uv_aux = cross([0 0 1], x_uv_aux)/norm(cross([0 0 1], x_uv_aux)); 
z_uv_aux = cross(x_uv_aux, y_uv_aux);

R_aux_I = [x_uv_aux, y_uv_aux', z_uv_aux'];

% Tranform coordinates from auxiliary to inertial reference frame
r_crater_aux = [x_crater_aux y_crater_aux z_crater_aux]';

r_crater_I = R_aux_I*r_crater_aux; % coordinate transformations
end
