function [pos_estimate,LS_e,A_mat,z_mat] = LS_estimate(meas_azimuth,meas_elevation, r_crater_LOS_nom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LS_ESTIMATE Least Square estimate of Spacecraft position
%  Copmutes a Least Square Estimate of Spacecraft position using angles by 
%  creating a pseudolinear measurement model 
%
% INPUTS
%   azimuth_meas     - (rad) measured azimuth angle in LOS frame
%   elevation_meas   - (rad) measured elevation angle in LOS frame
%   r_crater_LOS_nom - (m) known positions of craters in LOS frame
%
% OUTPUTS
%   pos_estimate - (m) the position estimae of the SC in LOS frame
%   LS_e         - (m) magnitude of position error
%   A_mat        - system matrix used for WLS estimation
%   z_mat        - measurement matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_craters = length(meas_azimuth); % num of detected craters

% Initialize matrices used in estimate
A = zeros(2*num_craters, 3);
z = zeros(2*num_craters, 1);

for i = 1:num_craters % stack matrices

    % Calc intermedoate vectors
    u_theta = [sin(meas_azimuth(i));-cos(meas_azimuth(i));0]; % Eq. 48

    u_psi = [-cos(meas_azimuth(i))*sin(meas_elevation(i)); ...
             -sin(meas_azimuth(i))*sin(meas_elevation(i)); ...
             cos(meas_elevation(i))]; % Eq. 52

    % Stacking into A matrix A_n [2x3]
    A(2*i-1:2*i, :) = [u_theta'; u_psi'];

    % Stacking into z vector z_n [2x1]
    p_n = r_crater_LOS_nom(:,i); % known position of considered crater 
    z(2*i-1:2*i) = [u_theta'*p_n; ...
                    u_psi'*p_n]; % measurement model

end

% Calc LS estimate of SC position
pos_estimate = (A'*A)\A'*z; % (m) in LOS frame
pos_estimate(2) = pos_estimate(2);

% Store absolute error
LS_e = norm(pos_estimate);

% Store system values
A_mat = A;
z_mat = z;

end