function [pos_estimate,BCWLS_e] = BCWLS_estimate(std_azimuth, std_elevation, meas_azimuth,meas_elevation,pos_WLS_est,A,W_block,do)
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
%   pos_WLS_est      - (m) WLS position estimate of SC in LOS frame
%   A                - system matrix from LS computation
%   W_block          - matrix of stacked weights
%   do               - (m) vector of distances from craters to LS estimate
%
% OUTPUTS
%   pos_estimate - (m) the position estimae of the SC in LOS frame
%   BCWLS_e      - (m) magnitude of position error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_craters = length(meas_azimuth); % num of detected craters

% Pre allocate m matrix
m_matrix = zeros(3,num_craters);

for i = 1:num_craters

    % Extract Wn and elements
    W_n = W_block(2*i-1:2*i, 2*i-1:2*i);
    a_n = W_n(1,1);
    b_n = W_n(1,2);
    c_n = W_n(2,2);

    % Form a vectors
    a_1n = [cos(meas_azimuth(i));sin(meas_azimuth(i));0];
    a_2n = [sin(meas_azimuth(i))*sin(meas_elevation(i));-cos(meas_azimuth(i))*sin(meas_elevation(i));0];
    a_3n = -[cos(meas_azimuth(i))*cos(meas_elevation(i));sin(meas_azimuth(i))*cos(meas_elevation(i));sin(meas_elevation(i))];

    % Form g vectors
    g_1n = -do(i)*cos(meas_elevation(i))*std_azimuth(i)^2*(a_n*a_1n+b_n*a_2n);
    g_2n = c_n*do(i)*a_3n*std_elevation(i)^2;

    % Calc m and save
    m_n = g_1n+g_2n;
    m_matrix(:,i) = m_n;
end

Atwn = sum(m_matrix,2);

gamma_gWLS = (A'*W_block*A)\Atwn;

pos_estimate = pos_WLS_est+gamma_gWLS;

BCWLS_e = norm(pos_estimate);
end