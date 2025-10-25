clear; close all; clc;

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CRATER_POS_CALC
%  - This script calculaes the position of a detected crater defined in
%  an inertial reference frame centered about the Moon as a function of the
%  position of the space craft and the incidence angle between the detcted
%  crater and spacecraft
%
% INPUTS
%  - distance: the altitude above the lunar surface (m)
%  - sc_incidence: incidence angle of spacecraft and crater (deg)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Parameters
radius_Moon = 1.7374e6; %(m)

% Inputs
distance = 760e3 + radius_Moon; % (m) from center of Moon
sc_incidence = 75*pi/180; % (rad) detected crater 


% Calculate position vector of spacecraft 
r_sc_I = [0.4811; 0.2580; 0.8379]*(distance); % (m)

% Calculate angle between r_sc_I and r_sc_I - r_c_I (alpha) using law of sines
d_sc = norm(r_sc_I); % magnitude
alpha = asin(radius_Moon*sin(pi-sc_incidence)/d_sc);

% Calculate angle between r_sc_I and r_c_I (gamma) 
beta = 180*pi/180 - sc_incidence; % angle between r_c_I and r_sc_I - r_c_I
gamma = pi - beta - alpha; % (rad)

% Define auxilary reference frame where gamma is the angle between the sc
% and crater from the center of the Moon
theta_circle = 0:0.01:2*pi; 
rand_int = randi(length(theta_circle)); % selects index of random angle theta
theta_rand = theta_circle(rand_int); % (rad) get random angle 
x_circle_aux = cos(gamma)*radius_Moon*ones(length(theta_circle), 1); % used to plot circle of possibilites
y_circle_aux = sin(gamma)*cos(theta_circle)*radius_Moon; 
z_circle_aux = sin(gamma)*sin(theta_circle)*radius_Moon;

x_crater_aux = cos(gamma)*radius_Moon; % distance along vector pointing to sc
y_crater_aux = sin(gamma)*cos(theta_rand)*radius_Moon; 
z_crater_aux = sin(gamma)*sin(theta_rand)*radius_Moon;

% Generate rotation matrix using unit vectors for auxiliary axes 
x_uv_aux = r_sc_I/d_sc;
y_uv_aux = cross([0 0 1], x_uv_aux)/norm(cross([0 0 1], x_uv_aux)); 
z_uv_aux = cross(x_uv_aux, y_uv_aux);

R_aux_I = [x_uv_aux, y_uv_aux', z_uv_aux'];

% Tranform coordinates from auxiliary to inertial reference frame
r_circle_aux = [x_circle_aux'; y_circle_aux; z_circle_aux]; % store positions in matrix
r_crater_aux = [x_crater_aux; y_crater_aux; z_crater_aux];

r_circle_I = R_aux_I*r_circle_aux; % coordinate transformations
r_crater_I = R_aux_I*r_crater_aux;

% Moon sphere coordinates
[xsphere, ysphere, zsphere] = sphere(30);
xsphere = radius_Moon * xsphere;
ysphere = radius_Moon * ysphere;
zsphere = radius_Moon * zsphere;

% plot in the auxiliary reference frame
figure;
axis equal
moon_sphere = surf(xsphere, ysphere, zsphere);
set(moon_sphere, 'FaceColor', [0.5 0.5 0.5]); 
hold on 
plot3(x_crater_aux, y_crater_aux, z_crater_aux, 'r*', 'MarkerSize',5, 'LineWidth',1.5)
plot3(x_circle_aux, y_circle_aux, z_circle_aux, 'm-', 'MarkerSize',5, 'LineWidth',2)
plot3(d_sc, 0, 0, 'b*', 'MarkerSize',5, 'LineWidth',2)
xlabel('x-axis (m)')
ylabel('y-axis (m)')
zlabel('z-axis (m)')
title('Auxiliary reference frame')
text(x_crater_aux, y_crater_aux, z_crater_aux, 'Crater', 'Color',[1 1 1])
text(d_sc, 0, 0, 'SC', 'Color',[1 1 1])

% plot in the inertial reference frame
figure;
axis equal
moon_sphere = surf(xsphere, ysphere, zsphere);
set(moon_sphere, 'FaceColor', [0.5 0.5 0.5]); 
hold on 
plot3(r_crater_I(1,:), r_crater_I(2,:), r_crater_I(3,:), 'yx', 'MarkerSize',5, 'LineWidth',1.5)
plot3(r_circle_I(1,:), r_circle_I(2,:), r_circle_I(3,:), 'g-', 'MarkerSize',5, 'LineWidth',1.5)
plot3(r_sc_I(1),r_sc_I(2), r_sc_I(3), 'b*', 'MarkerSize',5, 'LineWidth',2)
xlabel('x-axis (m)')
ylabel('y-axis (m)')
zlabel('z-axis (m)')
title('Inertial reference frame centered at the Moon')
text(r_crater_I(1), r_crater_I(2), r_crater_I(3), 'Crater', 'Color',[1 1 1])
text(r_sc_I(1), r_sc_I(2), r_sc_I(3), 'SC', 'Color',[1 1 1])
