close all; clear; clc;

altitude_sc = 35000; % current distance away from moon surface (need units)
moon_angle = 90; % incidence angle with moon surface (degrees)
radius_Moon = 1.7374e6; % m
distance_sc = altitude_sc+radius_Moon; % m
sc_bearing = 30; % deg


% Calculate position vector of spacecraft used to locate craters
r_sc_I = [0.4811; 0.2580; 0.8379]*(distance_sc); % m

% The following quantities were taken from discretized_results.mat for a
% crater with sc_inc = 30 deg and radius = 20
std_dev = 0.032129387741949; % rad
mean = 0.068789452907618; % rad

% Expand parameters for number of desired samples
N = 1000; % number of samples to be taken
std_dev_expanded = ones(N,1)*std_dev;
mean_expanded = ones(N,1)*mean;

[r_crater_I, r_crater_aux] = crater_pos(r_sc_I, sc_bearing, std_dev_expanded, mean_expanded);