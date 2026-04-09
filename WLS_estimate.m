function [pos_estimate,WLS_e,W_block,do] = WLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation, r_crater_LOS_nom, pos_LS_est,A,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WLS_ESTIMATE Weighted Least Square estimate of Spacecraft position
%  Copmutes a Weighted Least Square Estimate of Spacecraft position using 
%  angles by creating a pseudolinear measurement model 
%
% INPUTS
%   std_elevation    - (rad) analytical standard deviation of elevation
%   std_azimuth      - (rad) analytical standard deviation of azimuth
%   meas_azimuth     - (rad) measured azimuth angle in LOS frame
%   meas_elevation   - (rad) measured elevation angle in LOS frame
%   r_crater_LOS_nom - (m) known positions of craters in LOS frame
%   pos_LS_est       - (m) LS position estimate of SC in LOS frame
%   A                - system matrix from LS computation
%   z                - measurement matrix from LS computation
%
% OUTPUTS
%   pos_estimate - (m) the position estimae of the SC in LOS frame
%   WLS_e        - (m) magnitude of position error
%   W_block      - matrix of stacked weights
%   do           - (m) vector of distances from craters to LS estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_craters = length(meas_azimuth); % num of detected craters

% Calc distances of craters using LS estimate
do = vecnorm(pos_LS_est-r_crater_LOS_nom); % m

% Pre allocate spcae for weights Matrices
W_block = zeros(2*num_craters, 2*num_craters);

for i = 1:num_craters
    % Calc inter vectors using measured values
    uo_theta = [sin(meas_azimuth(i));-cos(meas_azimuth(i));0]; % Eq. 48

    uo_psi = [-cos(meas_azimuth(i))*sin(meas_elevation(i)); ...
             -sin(meas_azimuth(i))*sin(meas_elevation(i)); ...
             cos(meas_elevation(i))]; % Eq. 52
    
    % System matrix
    Ao = [uo_theta'; uo_psi']; 

    % Position error covariance
    R_p = zeros(3);

    % Measurement covariance
    R_theta = diag([std_azimuth(i)^2, std_elevation(i)^2]);

    % Temp matrix Do
    Do = diag([-do(i)*cos(meas_elevation(i)), do(i)]);

    % Compute covariance matrix R_eta
    R_eta_n = Do*R_theta*Do' + Ao*R_p*Ao'; % Eq. 55

    % Calc weight matrix
    W_n = inv(R_eta_n);

    W_block(2*i-1:2*i, 2*i-1:2*i) = W_n;
end

pos_estimate = (A'*W_block*A)\A'*W_block*z;
pos_estimate(2) = pos_estimate(2);

% Store absolute error
WLS_e = norm(pos_estimate);
end