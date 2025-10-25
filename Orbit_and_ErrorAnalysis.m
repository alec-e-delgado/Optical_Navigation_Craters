clc; close all; clear;

altitude = 35000; % current distance away from moon surface (need units)
moon_angle = 90; % incidence angle with moon surface (degrees)
radius_Moon = 1.7374e6; % m
distance_sc = altitude+radius_Moon; % m

% Calculate position vector of spacecraft used to locate craters
r_sc_I = [0.4811; 0.2580; 0.8379]*(distance_sc); % m

[angular_errors,sc_c_bearing,x_2, y_2, z_2] = angular_error_calc(altitude, moon_angle);

[r_craters_I, r_craters_aux] = crater_pos(r_sc_I, sc_c_bearing);

% --- Initialize figure for 3D error cone plot ---
clf; hold on; view(3);

% --- Plot each error vector from origin ---
i = 1;
while i <= length(x_2)
    plot3([0, x_2(i)], [0, y_2(i)], [0, z_2(i)], 'b');
    i = i + 1;
end

% --- Draw boundary circle of maximum deviation at x=1 ---
all_yz = [y_2; z_2];
r = max(abs(all_yz));
theta_circle = linspace(0, 2*pi, 100);
y = r .* cos(theta_circle);
z = r .* sin(theta_circle);
x = ones(size(theta_circle));
grid on;
plot3(x, y, z, 'b', 'LineWidth', 2);
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
xlim([-0.5,1.5]); ylim([-0.5,0.5]); zlim([-0.5,0.5]);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%% Propagating Orbit %%%%%%%%%%%%%%%%%%%%%%%

% --- Define basic orbital elements for a circular orbit ---
mu    = 4.9048695*10^12;     % gravitational parameter [m^3/s^2]
h     = 200*10^3;            % orbit altitude above surface [m]
a     = 1.74e6 + h;          % semi-major axis [m]
e     = 0;                   % eccentricity (0 → circular)
i     = deg2rad(30);         % inclination [rad]
LAN   = deg2rad(20);         % RA of ascending node [rad]
w     = deg2rad(0);          % argument of periapsis [rad]
M0    = deg2rad(10);         % initial mean anomaly [rad]

% --- Set up propagation timespan ---
t0      = 0;                         % start time [s]
t_final = 10 * 24 * 3600;            % 10 days → seconds
dt      = 100;                       % time step [s]
tspan   = [t0, t_final];             % for ODE45 solver

% --- Prepare time vector and compute mean anomaly array ---
t = t0:dt:t_final;                   % discrete time points [s]
N = numel(t);                        % number of steps
M = zeros(1, N);                     % preallocate mean anomaly
for k = 1:N
    delta_t = t(k) - t0;             % elapsed time [s]
    M(k)    = mod(M0 + sqrt(mu/a^3)*delta_t, 2*pi);  % Kepler's mean anomaly
end

% --- Solve Kepler's equation for eccentric anomaly E ---
E = M;                               % initial guess E₀ = M
for k = 1:N
    for iter = 1:10                  % Newton–Raphson iterations
        f_E       = E(k) - e*sin(E(k)) - M(k);
        f_prime_E = 1 - e*cos(E(k));
        E(k)      = E(k) - f_E / f_prime_E;  % update step
    end
end

% --- Compute true anomaly ν and radius r in orbital plane ---
nu = 2 * atan2(sqrt(1+e).*sin(E/2), sqrt(1-e).*cos(E/2));
rc = a * (1 - e*cos(E));             % orbital radius [m]

% --- Position & velocity in the orbital frame ---
ox = rc .* cos(nu);                  % x-coordinate
oy = rc .* sin(nu);                  % y-coordinate
oz = zeros(1, N);                    % orbital plane → z = 0

vx = -sqrt(mu/a) * sin(E);           
vy =  sqrt(mu/a) * sqrt(1 - e^2) .* cos(E);
vz = zeros(1, N);                    % orbital plane → vz = 0

orbital_pos = [ox; oy; oz];          
orbital_vel = [vx; vy; vz];

% --- Rotate into inertial frame via 3 Euler rotations ---
R1 = @(th)[cos(th) sin(th) 0; -sin(th) cos(th) 0; 0 0 1];  % about Z
R2 = @(th)[1 0 0; 0 cos(th) sin(th); 0 -sin(th) cos(th)];  % about X

R_total = R1(-LAN) * R2(-i) * R1(-w);
inertial_pos = R_total * orbital_pos;  
inertial_vel = R_total * orbital_vel;

% --- Convert to astronomical units (AU) and AU/day ---
AU = 1.49597870691e11;              % meters per AU
sec_per_day = 86400;               
r_AU   = inertial_pos / AU;        
v_AUpd = inertial_vel / (AU/sec_per_day);

% --- Plot the Keplerian orbit in 3D (start=green, end=red) ---
figure(2);
plot3(r_AU(1,:), r_AU(2,:), r_AU(3,:), 'b'); hold on;
plot3(r_AU(1,1), r_AU(2,1), r_AU(3,1), 'go', 'MarkerSize',8,'MarkerFaceColor','g');
plot3(r_AU(1,end), r_AU(2,end), r_AU(3,end), 'ro', 'MarkerSize',8,'MarkerFaceColor','r');
grid on; axis equal;


% Propagating orbit using ODE45 solver
phi0 = eye(6);                          % Initial state transition matrix (identity)

x0 = [inertial_pos(:,1);               % Initial position vector
      inertial_vel(:,1);               % Initial velocity vector
      reshape(phi0,36,1)];            % Flattened STM appended to state vector

% Set integration tolerances for accuracy
tol = odeset('RelTol',1e-8,'AbsTol',1e-9);

% Integrate equations of motion with STM propagation
[t_ode, X_ode] = ode45(@(t,x) eom_function(t, x, mu), tspan, x0, tol);

% Extract propagated state (position & velocity)
r_ode = X_ode(:,1:3)';                % Position solution
v_ode = X_ode(:,4:6)';                % Velocity solution

% Convert position to astronomical units for comparison
r_AU_ode = r_ode / AU;

% Reconstruct final state transition matrix from output
phi_final_flat = X_ode(end, 7:end);   % Flattened STM at final time
phi_final = reshape(phi_final_flat, 6, 6);  % STM in matrix form

% Interpolate ODE results onto Keplerian time grid for direct comparison
r_AU_ode_interp = zeros(3, length(t));
for i = 1:3
    r_AU_ode_interp(i,:) = interp1(t_ode, r_AU_ode(i,:), t, 'linear');  % Linear interpolation per axis
end

% Plot Keplerian vs. ODE45 trajectories
gcf;
plot3(r_AU(1,:), r_AU(2,:), r_AU(3,:), 'r', 'LineWidth', 1.5); hold on;
plot3(r_AU_ode(1,:), r_AU_ode(2,:), r_AU_ode(3,:), 'b--', 'LineWidth', 1.5);

% Mark start (green) and end (red) points on orbit
plot3(r_AU(1,1), r_AU(2,1), r_AU(3,1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot3(r_AU(1,end), r_AU(2,end), r_AU(3,end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');


% Title and legend
title('Comparison of Keplerian and ODE45 Orbital Trajectories');
legend('ODE45 Trajectory', 'Start Point', 'End Point', 'Keplerian Trajectory');
grid on; axis equal;

% Compute coordinate differences between methods
diff_x = r_AU_ode_interp(1,:) - r_AU(1,:);
diff_y = r_AU_ode_interp(2,:) - r_AU(2,:);
diff_z = r_AU_ode_interp(3,:) - r_AU(3,:);
time_days = t / 86400;                % Convert time to days for plotting

% Plot X-differences over time
figure;
plot(time_days, diff_x, 'b', 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('Difference in X'); grid on;

% Plot Y-differences over time
figure;
plot(time_days, diff_y, 'r', 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('Difference in Y'); grid on;

% Plot Z-differences over time
figure;
plot(time_days, diff_z, 'g', 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('Difference in Z'); grid on;

% Define landing cone half-angle and convert to radians
beta_deg = 5;
beta = deg2rad(beta_deg);

% Moon radius in meters
R_moon = 1737400;  

% Spacecraft inertial position [m]
x_sc = 3000e3;
y_sc = 1000e3;
z_sc = 500e3;

% Assemble spacecraft position vector
r_sc = [x_sc; y_sc; z_sc]; 

% Compute projection angles for crater rim
gamma = pi + asin((sin(beta) * norm(r_sc)) / R_moon);
theta_c = pi - gamma - beta;

% Generate random angles for crater points
rand_crater = 2 * pi * rand(1,30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute crater rim coordinates in spacecraft frame
x_crater_sc = R_moon * cos(theta_c) * ones(size(rand_crater));
y_crater_sc = R_moon * sin(theta_c) .* cos(rand_crater);
z_crater_sc = R_moon * sin(theta_c) .* sin(rand_crater);

% Build spacecraft local unit axes
x_sc_unit = r_sc / norm(r_sc);                                           % Radial direction
y_sc_unit = (cross([0 ;0 ;1], x_sc_unit)) / norm(cross([0 ;0 ;1], x_sc_unit)); % Perpendicular to z\ nz_sc_unit = cross(x_sc_unit, y_sc_unit);                               
z_sc_unit = cross(x_sc_unit, y_sc_unit);% Completes right-handed frame

% Form rotation matrix from spacecraft to inertial frame
R_sc2I = [x_sc_unit, y_sc_unit, z_sc_unit];

% Transform crater points to inertial coordinates
r_crater_I = R_sc2I * [x_crater_sc; y_crater_sc; z_crater_sc];

% Plot the Moon as a sphere
[xsphere, ysphere, zsphere] = sphere(30);
xsphere = R_moon * xsphere;
ysphere = R_moon * ysphere;
zsphere = R_moon * zsphere;
figure;
surf(xsphere, ysphere, zsphere)
hold on

% Plot spacecraft location
plot3(r_sc(1), r_sc(2), r_sc(3), 'ro');
hold on;

% Plot crater rim points
plot3(r_crater_I(1, :), r_crater_I(2, :), r_crater_I(3, :), 'b*');

axis equal;
grid on;

% Compute relative vector from spacecraft to first crater point
r_rel = r_crater_I(:,1) - r_sc;

% Build rotated inertial axes aligned with crater direction
x_r_unit = r_rel / norm(r_rel);
y_r_unit = cross([0; 0; 1], x_r_unit);
y_r_unit = y_r_unit / norm(y_r_unit);
z_r_unit = cross(x_r_unit, y_r_unit);

% Rotation matrix aligned with crater direction
R_I_rot = [x_r_unit, y_r_unit, z_r_unit];

% Example scaling vector (using provided variables)
r_2 = [x_2(1,1); y_2(1,1); y_2(1,1)];

% Scale rotated axes by r_2 values
R_2I = R_I_rot .* r_2;

% Plot landing cone rays based on angular error samples
for i = 1:length(angular_errors)
    r_cone_rel = [x_2(i); y_2(i); z_2(i)];                            % Cone direction in local frame
    r_cone_inertial = r_sc + (R_I_rot * r_cone_rel) * norm(r_rel);    % Transform to inertial
    plot3([r_sc(1) r_cone_inertial(1)], [r_sc(2) r_cone_inertial(2)], [r_sc(3) r_cone_inertial(3)], 'g-');
end

% Final plot settings
title('Spacecraft Position, Crater Rim, and Landing Cone Directions');
grid on;
axis equal;
view(3);
hold off;

% plot of craters 

% plot in the auxiliary reference frame
figure
axis equal
moon_sphere = surf(xsphere, ysphere, zsphere);
set(moon_sphere, 'FaceColor', [0.5 0.5 0.5]); 
hold on 
plot3(r_craters_aux(1,:), r_craters_aux(2,:), r_craters_aux(3,:), 'rx', 'MarkerSize',2, 'LineWidth',1.5)
plot3(distance_sc, 0, 0, 'b*', 'MarkerSize',5, 'LineWidth',2)
xlabel('x-axis (m)')
ylabel('y-axis (m)')
zlabel('z-axis (m)')
title('Auxiliary reference frame')
text(distance_sc, 0, 0, 'SC', 'Color',[1 1 1])

% plot in the inertial reference frame
figure;
axis equal
moon_sphere = surf(xsphere, ysphere, zsphere);
set(moon_sphere, 'FaceColor', [0.5 0.5 0.5]); 
hold on 
plot3(r_craters_I(1,:), r_craters_I(2,:), r_craters_I(3,:), 'yx', 'MarkerSize',2, 'LineWidth',1.5)
plot3(r_sc_I(1),r_sc_I(2), r_sc_I(3), 'b*', 'MarkerSize',5, 'LineWidth',2)
xlabel('x-axis (m)')
ylabel('y-axis (m)')
zlabel('z-axis (m)')
title('Inertial reference frame centered at the Moon')
text(r_sc_I(1), r_sc_I(2), r_sc_I(3), 'SC', 'Color',[1 1 1])


% Equations of motion and STM propagation for two-body problem
function dxdt = eom_function(t, x, mu)

    % Unpack state vector into position (r), velocity (v), and STM (phi)
    r = x(1:3);
    v = x(4:6);
    phi = reshape(x(7:42), 6, 6);

    % Compute norm of position and gravitational acceleration
    r_norm = norm(r);
    a = -mu / r_norm^3 * r;

    % Extract position components and compute r^3 and r^5
    x_pos = r(1); y_pos = r(2); z_pos = r(3);
    r3 = r_norm^3;
    r5 = r_norm^5;

    % Build matrix of partial derivatives of acceleration w.r.t. position
    dadr = mu * [ (3*x_pos^2)/r5 - 1/r3, (3*x_pos*y_pos)/r5,       (3*x_pos*z_pos)/r5;
                  (3*x_pos*y_pos)/r5,    (3*y_pos^2)/r5 - 1/r3,    (3*y_pos*z_pos)/r5;
                  (3*x_pos*z_pos)/r5,    (3*y_pos*z_pos)/r5,       (3*z_pos^2)/r5 - 1/r3  ];

    % Assemble system matrix A for variational equations
    A = [zeros(3), eye(3); dadr, zeros(3)];

    % Compute time derivative of STM: dphi/dt = A * phi
    dphi = A * phi;

    % Pack derivatives: velocity, acceleration, and flattened STM
    dxdt = [v; a; reshape(dphi, 36, 1)];
end



