close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EXTENDED KALMAN FILTER SCRIPT (EKF) SCRIPT
% This script implements an EKF over a lunar orbit for position estimation
%
% DESCRIPTION
%  - Spacecraft (sc) is in a lunar orbit and the LS main script is used to
%    initialize an EKF, the EKF is used to refine position estimations over
%    the orbit of the sc
%
% EXPERIMENTS
%  - Circular orbit: circular orbit about the moon at constant solar phase
%    angle
%
% INPUTS
%  - Initial conditions for orbits
%
% OUTPUTS
%  - Plotted simulation results
%
% Author: Alec Delgado
% Date: 2026-03-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input

% Time
t0 = 0.1; 
tf = 300; 
tspan_meas = t0:0.1:tf; % for measurement generation
N = length(tspan_meas); % number of measurements

% Circular orbit experiment


% Initial EKF estimate
xhat0 = 0.8 * x0_true;

% Initial covariance
P0 = 10e10 * eye(2);

% ODE options
options = odeset('RelTol',1e-12,'AbsTol',1e-9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Task 6.1

% Simulate truth and measurement
[x_true, z] = simulate_truth_and_measurements(tspan_meas, xhat0, p_true, sigma_v, options);

% Process Covariance
q62 = 1e-6;
Q61 = q62 * eye(2);

p_est_61 = p_true;

[xhat_61, P_61] = run_cd_ekf(tspan_meas, z, xhat0, P0, p_est_61, Q61, R, H, options);

err_61 = xhat_61 - x_true;
sig_61 = zeros(2,N);
for k = 1:N
    sig_61(:,k) = sqrt(diag(P_61(:,:,k)));
end

%% Task 6.2

% Ignore p4-p6 in estimate
p_est_62 = p_true;
p_est_62(4:6) = 0;   

% Process Covariance
q62 = 0.05;
Q62 = q62 * eye(2);

[xhat_62, P_62] = run_cd_ekf(tspan_meas, z, xhat0, P0, p_est_62, Q62, R, H, options);

err_62 = xhat_62 - x_true;
sig_62 = zeros(2,N);
for k = 1:N
    sig_62(:,k) = sqrt(diag(P_62(:,:,k)));
end


%% 6.1 plots
figure % Truth and measurements overlaid
subplot(2,1,1)
hold on
plot(tspan_meas, x_true(1,:), 'k', 'LineWidth', 1.5)
plot(tspan_meas, xhat_61(1,:), 'r--', 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('x')
legend('True','Estimated','Location','best')
box on;

subplot(2,1,2)
hold on
plot(tspan_meas, x_true(2,:), 'k', 'LineWidth', 1.5)
plot(tspan_meas, xhat_61(2,:), 'r--', 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('$\dot{x}$','Interpreter','latex')
legend('True','Estimated','Location','best')
box on;

figure % Errors and associated bounds
subplot(2,1,1);
hold on
plot(tspan_meas(2:end), err_61(1,2:end), 'b', 'LineWidth', 1.1)
plot(tspan_meas(2:end),  3*sig_61(1,2:end), 'r--', 'LineWidth', 1.1)
plot(tspan_meas(2:end), -3*sig_61(1,2:end), 'r--', 'LineWidth', 1.1)
grid on
xlabel('Time (s)')
ylabel('x error')
legend('error','+3\sigma','-3\sigma','Location','best')
box on;

subplot(2,1,2)
hold on
plot(tspan_meas(3:end), err_61(2,3:end), 'b', 'LineWidth', 1.1)
plot(tspan_meas(3:end),  3*sig_61(2,3:end), 'r--', 'LineWidth', 1.1)
plot(tspan_meas(3:end), -3*sig_61(2,3:end), 'r--', 'LineWidth', 1.1)
grid on
xlabel('Time (s)')
ylabel('$\dot{x}$ error','Interpreter','latex')
legend('error','+3\sigma','-3\sigma','Location','best')
box on;

%% 6.2 plots

figure % Truth and measurements overlaid
subplot(2,1,1)
hold on
plot(tspan_meas, x_true(1,:), 'k', 'LineWidth', 1.5)
plot(tspan_meas, xhat_62(1,:), 'r--', 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('x')
legend('True','Estimated','Location','best')
box on;

subplot(2,1,2)
hold on
plot(tspan_meas, x_true(2,:), 'k', 'LineWidth', 1.5)
plot(tspan_meas, xhat_62(2,:), 'r--', 'LineWidth', 1.2)
grid on
xlabel('Time (s)')
ylabel('$\dot{x}$','Interpreter','latex')
legend('True','Estimated','Location','best')
box on;

figure % Errors and associated bounds
subplot(2,1,1);
hold on
plot(tspan_meas(2:end), err_62(1,2:end), 'b', 'LineWidth', 1.1)
plot(tspan_meas(2:end),  3*sig_62(1,2:end), 'r--', 'LineWidth', 1.1)
plot(tspan_meas(2:end), -3*sig_62(1,2:end), 'r--', 'LineWidth', 1.1)
grid on
xlabel('Time (s)')
ylabel('x error')
legend('error','+3\sigma','-3\sigma','Location','best')
box on;

subplot(2,1,2)
hold on
plot(tspan_meas(3:end), err_62(2,3:end), 'b', 'LineWidth', 1.1)
plot(tspan_meas(3:end),  3*sig_62(2,3:end), 'r--', 'LineWidth', 1.1)
plot(tspan_meas(3:end), -3*sig_62(2,3:end), 'r--', 'LineWidth', 1.1)
grid on
xlabel('Time (s)')
ylabel('$\dot{x}$ error','Interpreter','latex')
legend('error','+3\sigma','-3\sigma','Location','best')
box on;

%% Helper functions

function [x_true, z] = simulate_truth_and_measurements(t_meas, x0_true, p_true, sigma_v, opts)
    % Simulate truth over full interval
    [t_true, X_true] = ode45(@(t,s) true_dynamics(t, s, p_true), [t_meas(1), t_meas(end)], x0_true, opts);

    % Sample at measurement times
    x_true = interp1(t_true, X_true, t_meas, 'pchip').';  % 2 x N

    % Noisy position measurements
    z = x_true(1,:) + normrnd(0,sigma_v,1,length(t_meas));
end

function ds = true_dynamics(t, s, p)
    % True nonlinear system
    % x(1)=position, x(2)=velocity
    p1 = p(1); p2 = p(2); p3 = p(3);
    p4 = p(4); p5 = p(5); p6 = p(6);

    ds = zeros(2,1);
    ds(1) = s(2);
    ds(2) = -p1*s(2) - p2*s(1) - p3*s(1)^3 - p4*sin(p5*t + p6);
end

function dX = ekf_aug_ode(t, X, p, Q)
    % Augmented ODE for EKF propagation:
    xhat = X(1:2);
    P = reshape(X(3:6), 2, 2);

    % Jacobian df/dX
    x = X(1);

    p1 = p(1);
    p2 = p(2);
    p3 = p(3);

    A = [0, 1;
        -(p2 + 3*p3*x^2), -p1];

    dxhat = true_dynamics(t, xhat, p);
    Pdot = A*P + P*A' + Q;

    dX = [dxhat; Pdot(:)];
end

function [xhat_hist, P_hist] = run_cd_ekf(t_meas, z, xhat0, P0, p_est, Q, R, H, opts)
    N = length(t_meas);

    xhat_hist = zeros(2,N);
    P_hist = zeros(2,2,N);

    xhat_plus = xhat0;
    P_plus = P0;

    xhat_hist(:,1) = xhat_plus;
    P_hist(:,:,1) = P_plus;

    for k = 2:N
        tk_1 = t_meas(k-1);
        tk = t_meas(k);

        % Continuous propagation of state and covariance
        X0 = [xhat_plus; P_plus(:)];
        [~, Xprop] = ode45(@(t,X) ekf_aug_ode(t, X, p_est, Q), [tk_1 tk], X0, opts);

        Xminus = Xprop(end,:)';
        xhat_minus = Xminus(1:2);
        P_minus = reshape(Xminus(3:6), 2, 2);

        % enforce symmetry
        P_minus = 0.5*(P_minus + P_minus.');

        % measurement update
        S = H*P_minus*H' + R;
        K = (P_minus*H') / S;

        innov = z(k) - H*xhat_minus;
        xhat_plus = xhat_minus + K*innov;

        I2 = eye(2);
        P_plus = (I2 - K*H)*P_minus;   % Joseph form
        P_plus = 0.5*(P_plus + P_plus.');

        xhat_hist(:,k) = xhat_plus;
        P_hist(:,:,k) = P_plus;
    end
end