clc;
clearvars;
close all;

num_drums =4; % select number of drums 8,4,2,1

dt = 0.1;
T = 6000;
n0 = 1;
time = 0:dt:T;
nt = length(time);
ref = zeros(nt, 1);



function [Rho_d0, Reactivity_per_degree, u0] = number_drums(num_drums)
    switch num_drums
        case 8
            Rho_d0 = -0.033085599;
            Reactivity_per_degree = 26.11e-5;
            u0 = 77.56; % 8drums 
        case 4
            Rho_d0 = -0.033085599+0.013980296;
            Reactivity_per_degree = 16.11e-5; 
            u0 = 108.5; % 4 drums 
        case 2
            Rho_d0 = -0.033085599 + 0.0074; 
            Reactivity_per_degree = 7.33e-5;
            u0 = 165.5; % 2 drums
        case 1
            Rho_d0 = -0.033085599 + 0.0071 + 0.0082; 
            Reactivity_per_degree = 2.77e-5;
            u0 = 170;% 1 drums
        otherwise
            Rho_d0 = -0.033085599;
            Reactivity_per_degree = 26.11e-5;
            u0 = 77.56; % 8drums 
    end
end



time_point=[ 0 20 30 50 60 80 90 110 130 200]*30;

if num_drums == 8
        pow=[ 1 1 0.5 0.5 1 1 0.5 0.5 1 1]; 
    elseif num_drums == 4
        pow=[0.3 0.3 1 1 0.6 0.6 0.8 0.8 1 1];
    elseif num_drums == 2
        pow=[0.7 0.7 0.4 0.4 0.4 0.8 0.8 0.8 1 1]; 
    elseif num_drums == 1
        pow=[0.9 0.9 0.7 0.7 0.5 0.5 0.7 0.7 1 1]; 
    else
        error('Invalid number of drums specified.');
end
 

ref_old = pow(1);
for it = 1:nt 
    if (it > 1)
        time(it) = time(it - 1) + dt;
        ref_old = ref(it - 1);  
    end
    ref(it) = ref_old;
    i1 = 1;
    i2 = 1;
    for ii = 1:length(time_point) - 1
        if (time_point(ii) <= time(it) && time(it) <= time_point(ii + 1))
            i1 = ii;
            i2 = ii + 1;
            frac1 = (time_point(ii + 1) - time(it)) / (time_point(ii + 1) - time_point(ii));
            frac2 = 1.0 - frac1;
            ref(it) = frac1 * pow(i1) + frac2 * pow(i2);
            break
        end
    end
end

ref(:) = ref(:);
yref = (ref .* n0);

Kp = 2;
Ki = 5;
Kaw = 0.3;
Kd = 0.001;
T_C = 0.2;
max = 180;
min = 0;
max_rate = 0.5;

[Rho_d0, Reactivity_per_degree, u0] = number_drums(num_drums);

Kp = Kp * 26.11e-5 / Reactivity_per_degree;
Ki = Ki * 26.11e-5 / Reactivity_per_degree;
Kaw = Kaw * 26.11e-5 / Reactivity_per_degree;
Kd = Kd * 26.11e-5 / Reactivity_per_degree;
max_rate = max_rate * 26.11e-5 / Reactivity_per_degree;

%% initial condition
Sig_x = 2.65e-22;
yi = 0.061;
yx = 0.002;
lamda_x = 2.09e-5;
lamda_I = 2.87e-5;
Sum_f = 0.3358;
G = 3.2e-11;
V = 400 * 200;
P_0 = 22e6;
Pi = P_0 / (G * Sum_f * V);
Xe0 = (yi + yx) * Sum_f * Pi / (lamda_x + Sig_x * Pi);
I0 = yi * Sum_f * Pi / lamda_I;




% Select initial condition based on the number of drums
if num_drums == 8
        x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 900 898 883];
    elseif num_drums == 4
        x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 875 873.5 870];
    elseif num_drums == 2
        x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 890 888 877];
    elseif num_drums == 1
        x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 897 895 881];
    else
        error('Invalid number of drums specified.');
end




pidController(0, 0, 0, Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 1, u0);

%% Simulation
x(1, :) = x0;
u(1) = u0;
rho = zeros(nt, 1);
for i = 1:nt - 1
    [dx, rho(i)] = reactorDAE(0, x(i,:), u(i), Rho_d0, Reactivity_per_degree,Xe0, I0, Pi,num_drums);
    x(i + 1, :) = x(i, :) + dx' * dt;
    u(i + 1) = pidController(time(i), x(i + 1, 1), ref(i + 1), Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, 0, 0);
end
%[~, rho(nt)] = reactorDAE(0, x(nt, :), u(nt), Rho_d0, Reactivity_per_degree,Xe0, I0, Pi); 

 t = time;
% du = [0 diff(u) / dt];

du = zeros(size(u));  
du(1:2) = 0.1;        
du(3:end) = diff(u(2:end)) / dt; 
control_effort_withoutFilter=norm(du);
fprintf('control_effort_withoutFilter: %f\n', control_effort_withoutFilter);

window_size = 2000;  % Adjust this window size to control the degree of smoothing
du = movmean(du, window_size);  % Smooth du with a moving average filter because of the method of differention 


%% Interpolate ref in case t and time are different
ref_interpolated = interp1(time, ref, t, 'linear');
r = ref_interpolated;
y = x(:, 1);

%% Error analysis and display
% Downsampling
downsample_factor = 10; % Adjust this factor based on your dataset
t_ds = t(1:downsample_factor:end);
r_ds = r(1:downsample_factor:end);
y_ds = y(1:downsample_factor:end);

% Calculate the error
e_ds = r_ds - y_ds';

% Mean Absolute Error (MAE)
MAE = mean(abs(e_ds));

% Integral Absolute Error (IAE)
IAE = trapz(t_ds, abs(e_ds));

% Integral Time Absolute Error (ITAE)
ITAE = trapz(t_ds, t_ds .* abs(e_ds));

% Integral Square Error (ISE)
ISE = trapz(t_ds, e_ds.^2);

% Integral Time Square Error (ITSE)
ITSE = trapz(t_ds, t_ds .* e_ds.^2);

% Display results once
fprintf('Mean Absolute Error (MAE): %.4f\n', MAE);
fprintf('Integral Absolute Error (IAE): %.4f\n', IAE);
fprintf('Integral Time Absolute Error (ITAE): %.4f\n', ITAE);
fprintf('Integral Square Error (ISE): %.4f\n', ISE);
fprintf('Integral Time Square Error (ITSE): %.4f\n', ITSE);


control_effort = norm(du);
fprintf('Control Effort: %f\n', control_effort);

%% Plot results


figure();

% Adjust the font size and font weight to bold
font_size = 12;
font_weight = 'bold';
t=t/60;
subplot(4,2,[1,2]);
plot(t, x(:, 1) * 100, 'LineWidth', 2);
hold on;
plot(t, ref_interpolated * 100, '--', 'LineWidth', 2);
hold off;
grid on;
xlabel('Time (min)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Relative Power (%)', 'FontSize', font_size, 'FontWeight', font_weight);
title('PID microreactor control system core power simulation with 8 control drums', 'FontSize', font_size, 'FontWeight', font_weight);
ylim([0, 120]);
xlim([0,100])
legend({'Actual power', 'Desired power'}, 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); 

subplot(4,2,3);
plot(t, u, 'LineWidth', 2);
grid on;
xlabel('Time (min)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Rotation (deg)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Drum Rotation', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); 
xlim([0, 100]);

% switch num_drums
%     case 8
%         ylim([65, 85]);
%     case 4
%         ylim([110, 135]);
%     case 2
%         ylim([120, 185]);
%     case 1
%         ylim([120, 190]);
%     otherwise
%         ylim([65, 85]);
% end

subplot(4,2,4);
plot(t, du, 'LineWidth', 2);
grid on;
xlabel('Time (min)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('speed (deg/s)', 'FontSize', font_size, 'FontWeight', font_weight);
title('conrol drum speed', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold
xlim([0, 100]);

subplot(4,2,5);
plot(t, rho*1e5, 'LineWidth', 2);
%plot(t(500:end), rho(500:end)*1e5, 'LineWidth', 2);
grid on;
xlabel('Time (min)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Reactivity (PCM)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Core Reactivity', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold
xlim([0, 100]);

subplot(4,2,6);
plot(t, x(:, 10), 'LineWidth', 2);
%plot(t(700:end), x(700:end, 10), 'LineWidth', 2);
grid on;
xlabel('Time (min)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Temperature (K)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Temperature of fuel', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold
xlim([0, 100]);

subplot(4,2,7);
plot(t, x(:, 11), 'LineWidth', 2);
grid on;
xlabel('Time (min)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Temperature (K)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Temperature of moderator', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold
xlim([0, 100]);

subplot(4,2,8);
plot(t, x(:,12), 'LineWidth', 2);
grid on;
xlabel('Time (min)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Temperature (K)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Temperature of coolant', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold
xlim([0, 100]);

% Adjust the layout so labels don't overlap
set(gcf, 'Position', [5, 5, 2000, 4500]);

saveas(gcf, 'PID_power_simulation.png');

function [dx, rho] = reactorDAE(t, x, u, Rho_d0, Reactivity_per_degree, Xe0, I0, Pi,num_drums)
    %% Parameters
    Sig_x   = 2.65e-22;
    yi      = 0.061;
    yx      = 0.002;
    lamda_x = 2.09e-5;
    lamda_I = 2.87e-5;
    Sum_f   = 0.3358;
    
    l       = 1.68e-3;
    beta    = 0.0048;
    beta_1  = 1.42481E-04; 
    beta_2  = 9.24281E-04;
    beta_3  = 7.79956E-04;
    beta_4  = 2.06583E-03;
    beta_5  = 6.71175E-04;
    beta_6  = 2.17806E-04;
    Lamda_1 = 1.272E-02;
    Lamda_2 = 3.174E-02;
    Lamda_3 = 1.160E-01;
    Lamda_4 = 3.110E-01;
    Lamda_5 = 1.400E+00;
    Lamda_6 = 3.870E+00;

    cp_f    = 977;
    cp_m    = 1697;
    cp_c    = 5188.6;
    M_f     = 2002;
    M_m     = 11573;
    M_c     = 500;
    mu_f    = M_f * cp_f;
    mu_m    = M_m * cp_m;
    mu_c    = M_c * cp_c;
    f_f     = 0.96;
    P_0     = 22e6; 
    Tf0     = 1105;
    Tm0     = 1087;
    T_in    = 864;
    T_out   = 1106;
    Tc0     = (T_in + T_out) / 2;
     
    K_fm    = f_f * P_0 / (Tf0 - Tm0);
    K_mc    = P_0 / (Tm0 - Tc0);
    M_dot   = 1.75E+01;
    alpha_f = -2.875e-5;
    alpha_m = -3.696e-5;
    alpha_c = 0.0;

    %% Declaration of state variables, x(i), where i = 1 to 14
    n_r     = x(1); 
    Cr1     = x(2); 
    Cr2     = x(3); 
    Cr3     = x(4); 
    Cr4     = x(5);
    Cr5     = x(6);
    Cr6     = x(7);
    X       = x(8); 
    I       = x(9); 
    Tf      = x(10); 
    Tm      = x(11); 
    Tc      = x(12); 

    Rho_d1 = Rho_d0 + u * Reactivity_per_degree;
%900.42 898.28 888.261
    %% ODEs
    dx      = zeros(12, 1);
    %rho     = Rho_d1 + alpha_f * (Tf - Tf0) + alpha_c * (Tc - Tc0) + alpha_m * (Tm - Tm0) - Sig_x * (X - Xe0) / Sum_f; % use for 8,2,1 drums
    %rho     = Rho_d1 + alpha_f * (Tf - 900.42) + alpha_c * (Tc - 888.261) + alpha_m * (Tm - 898.261) - Sig_x * (X - Xe0) / Sum_f; %4drums

    if num_drums == 8 || num_drums == 2 || num_drums == 1
        rho = Rho_d1 + alpha_f * (Tf - Tf0) + alpha_c * (Tc - Tc0) + alpha_m * (Tm - Tm0) - Sig_x * (X - Xe0) / Sum_f;
    elseif num_drums == 4
        rho = Rho_d1 + alpha_f * (Tf - 900.42) + alpha_c * (Tc - 888.261) + alpha_m * (Tm - 898.261) - Sig_x * (X - Xe0) / Sum_f;
    else
        error('Invalid number of drums specified.');
    end
    %% Kinetics equations with six-delayed neutron groups
    dx(1)   = (rho - beta) / l * n_r + beta_1 / l * Cr1 + beta_2 / l * Cr2 + beta_3 / l * Cr3 + beta_4 / l * Cr4 + beta_5 / l * Cr5 + beta_6 / l * Cr6;
    dx(2)   = Lamda_1 * n_r - Lamda_1 * Cr1;
    dx(3)   = Lamda_2 * n_r - Lamda_2 * Cr2;
    dx(4)   = Lamda_3 * n_r - Lamda_3 * Cr3;
    dx(5)   = Lamda_4 * n_r - Lamda_4 * Cr4;
    dx(6)   = Lamda_5 * n_r - Lamda_5 * Cr5;
    dx(7)   = Lamda_6 * n_r - Lamda_6 * Cr6;
    
    %% Xenon and Iodine dynamics
    dx(8)   = yx * Sum_f * Pi + lamda_I * I - Sig_x * X * Pi - lamda_x * X;
    dx(9)   = yi * Sum_f * Pi - lamda_I * I;
    
    %% Thermal–hydraulics model of the reactor core
    dx(10)  = f_f * P_0 / mu_f * n_r - K_fm / mu_f * (Tf - Tc);
    dx(11)  = (1 - f_f) * P_0 / mu_m * n_r + (K_fm * (Tf - Tm) - K_mc * (Tm - Tc)) / mu_m;
    dx(12)  = K_mc * (Tm - Tc) / mu_c - 2 * M_dot * cp_c * (Tc - T_in) / mu_c;

    
end

%% PID controller
function u = pidController(t, measurement, setpoint, Kp, Ki, Kd, Kaw, T_C, max, min, max_rate, init, u0)

    persistent integral err_prev deriv_prev t_prev command_sat_prev command_prev;

    % Initialise static variables    
    if (isempty(integral) || init == 1)
        integral         = u0;
        err_prev         = 0;
        t_prev           = 0;
        deriv_prev       = 0;
        command_prev     = u0;
        command_sat_prev = u0;
    end
    
    % Calculate error
    err = setpoint - measurement;

    % Calculate period
    T = t - t_prev;

    % Store previous time
    t_prev = t;
    
    % Update integral term with anti-windup
    integral = integral + Ki * err * T + Kaw * (command_sat_prev - command_prev) * T;
    
    % Calculate filtered derivative
    deriv_filt = (err - err_prev + T_C * deriv_prev) / (T + T_C);

    % Store previous error and previous derivative
    err_prev = err;
    deriv_prev = deriv_filt;
    
    % Calculate command using PID equation
    command = Kp * err + integral + Kd * deriv_filt;
    
    % Store previous command
    command_prev = command;
    
    % Saturate command
    if command > max
        command_sat = max;
    elseif command < min
        command_sat = min;
    else
        command_sat = command;
    end
    
    % Apply rate limiter
    if command_sat > command_sat_prev + max_rate * T
        command_sat = command_sat_prev + max_rate * T;
    elseif command_sat < command_sat_prev - max_rate * T
        command_sat = command_sat_prev - max_rate * T;
    end
    
    % Store previous saturated command
    command_sat_prev = command_sat;
    
    % Output the saturated command
    u = command_sat;
end
