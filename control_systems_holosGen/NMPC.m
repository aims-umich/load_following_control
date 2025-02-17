clear
close all
clc

% Define constants and initial parameters
nx = 12;
ny = 1;
nu = 1;
dt = 0.1;  
Ts = 0.01; 

% Create NL-MPC object
nlobj = nlmpc(nx, ny, nu);
nlobj.Ts = Ts;

nlobj.PredictionHorizon = 10; 
nlobj.ControlHorizon = 3;    

% nlobj.PredictionHorizon = 5; 
% nlobj.ControlHorizon = 2;   

% Model setup
nlobj.Model.StateFcn = @reactorDT0;
nlobj.Model.IsContinuousTime = false;
nlobj.Model.NumberOfParameters = 1;

%define number of drums 
num_drums = 8;  % Define number of drums and define in the discrete function 

nlobj.Model.OutputFcn = @(x, u, params) reactorOutputFcn(x, u, params, num_drums);
nlobj.Jacobian.OutputFcn = @(x, u, Ts) [1 0 0 0 0 0 0 0 0 0 0 0];


% nlobj.Weights.OutputVariables = 1e7; 
% nlobj.Weights.ManipulatedVariables = 0.1;
% nlobj.Weights.ManipulatedVariablesRate = 1;

nlobj.Weights.OutputVariables = 1e7; 
nlobj.Weights.ManipulatedVariables = 0.1;
nlobj.Weights.ManipulatedVariablesRate = 1;



% Manipulated variable constraints
nlobj.MV.Min = 0;  
nlobj.MV.Max = 180; 
nlobj.MV.RateMin = -0.5 * dt;
nlobj.MV.RateMax = 0.5 * dt;
% nlobj.MV.RateMin = -1 * dt; %this works with 1 drum when power changing is high
% nlobj.MV.RateMax = 1 * dt;

% Simulation parameters
%T = 1000; % Reduced total simulation time
T = 6000;
time = 0:dt:T;
nt = length(time);
ref = zeros(nt, 1);

% Function to set parameters based on number of drums
[Rho_d0, Reactivity_per_degree, u0] = setParameters(num_drums);

% Time-varying reference for simulation

%time_point=[0 20 30 50 60 80 90 110 130 150];
%time_point=[0 10 40 70 80 100 110 130 140 150];
%time_point=[0 30 40 50 60 70 90 100 140 150];
time_point=[0 20 30 50 60 80 90 110 130 200]*30;
%pow=[0.8 0.8 0.4 0.4 0.8 0.8 0.4 0.4 0.8 0.8];
pow=[1 1 0.5 0.5 1 1 0.5 0.5 1 1]; % use for 8 drums 
%pow=[0.3 0.3 1 1 0.6 0.6 0.8 0.8 1 1];% use for 4 drums
%pow=[0.7 0.7 0.4 0.4 0.4 0.8 0.8 0.8 1 1]; %use for 2 drums
%pow=[0.9 0.9 0.7 0.7 0.5 0.5 0.7 0.7 1 1]; % use for 1 drum 
%pow=[1 1 1 1 1 1 1 1 1 1];
%pow=[1 1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];emergncy or power reduction


ref_old = pow(1);

% Generate reference trajectory
for it = 1:nt
    if(it > 1)
        time(it) = time(it - 1) + dt;
        ref_old = ref(it - 1);
    end
    ref(it) = ref_old;
    for ii = 1:length(time_point) - 1
        if(time_point(ii) <= time(it) && time(it) <= time_point(ii + 1))
            frac1 = (time_point(ii + 1) - time(it)) / (time_point(ii + 1) - time_point(ii));
            frac2 = 1.0 - frac1;
            ref(it) = frac1 * pow(ii) + frac2 * pow(ii + 1);
            break;
        end
    end
end

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


%x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 900.42 898.28 888.261]';

x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 900 898 883]';%8 drums
%x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 875 873.5 870]'; % 4 drums 
%x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 890 888 877]';%2 drums
%x0 = [pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) pow(1) I0 Xe0 897 895 881]';%1 drum

% Extended Kalman filter
EKF = extendedKalmanFilter(@reactorStateFcn, @reactorMeasurementFcn);
x = x0;
y = x0(1);
EKF.State = x0;
mv = u0;

% NL-MPC options
nloptions = nlmpcmoveopt;
nloptions.Parameters = {Ts};
nloptions.X0 = x0';
nloptions.MV0 = mv;

% Initialize simulation variables
Duration = 6000; 
hbar = waitbar(0, 'Simulation Progress');
xHistory = zeros(12, Duration / dt + 1);
xkHistory = zeros(12, Duration / dt + 1);
MV = zeros(1, Duration / dt + 1);

xHistory(:, 1) = x0;
xkHistory(:, 1) = x0;
MV(1) = u0;
rho = zeros(1, nt);
% Main simulation loop
for ct = 1:(Duration / dt)
    xk = correct(EKF, y);
    [mv, nloptions, info] = nlmpcmove(nlobj, xk, mv, ref(ct), [], nloptions);
    xk = predict(EKF, [mv; Ts]);
    x = reactorDT0(x, mv, Ts);
    y = x(1);
    [~, rho(ct)] = reactorCT0(x, mv, Rho_d0, Reactivity_per_degree);

    xHistory(:, ct) = x;
    xkHistory(:, ct) = xk;
    MV(ct) = mv;
    waitbar(ct * dt / Duration, hbar);
end
close(hbar)

du = [0 diff(MV) / dt];
x_final = xHistory(:, end);  
u_final = MV(end);  

% du = zeros(size(MV));  
% du(1:2) = 0.1;        
% du(3:end) = diff(MV(2:end)) / dt; 
% window_size = 2000;  
% du = movmean(du, window_size);   




[~, rho(nt)] = reactorCT0(xHistory(:, nt), MV(nt), Rho_d0, Reactivity_per_degree);

control_effort = norm(du);



fprintf('Control Effort: %f\n', control_effort);
%% Error analysis and display
error = xHistory(1, :) - ref';
t = 0:dt:Duration;
% Mean Absolute Error (MAE)
MAE = mean(abs(error));

% Integral Absolute Error (IAE)
IAE = trapz(t, abs(error));

% Integral Time Absolute Error (ITAE)
ITAE = trapz(t, t .* abs(error));

% Integral Square Error (ISE)
ISE = trapz(t, error.^2);

% Integral Time Square Error (ITSE)
ITSE = trapz(t, t .* error.^2);

% Display results once
fprintf('Mean Absolute Error (MAE): %.4f\n', MAE);
fprintf('Integral Absolute Error (IAE): %.4f\n', IAE);
fprintf('Integral Time Absolute Error (ITAE): %.4f\n', ITAE);
fprintf('Integral Square Error (ISE): %.4f\n', ISE);
fprintf('Integral Time Square Error (ITSE): %.4f\n', ITSE);


%% Plotting results

t = 0:dt:Duration;

figure();

% Adjust the font size and font weight to bold
font_size = 12;
font_weight = 'bold';
t=t/60;
subplot(4,2,[1,2]);
plot(t(1:end-1), xHistory(1, 1:end-1) * 100, 'LineWidth', 2)
hold on;
plot(t(1:end-1), ref(1:end-1) * 100, '--', 'LineWidth', 2);
hold off;
grid on;
xlabel('Time (mins)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Relative Power (%)', 'FontSize', font_size, 'FontWeight', font_weight);
title('NMPC microreactor control system core power simulation with 4 control drums', 'FontSize', font_size, 'FontWeight', font_weight);
ylim([0, 200]);
legend({'Actual power', 'Desired power'}, 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); 

subplot(4,2,3);
% plot(time, u, 'LineWidth', 2);
plot(t(1:end-1), MV(1:end-1), 'LineWidth', 2);
grid on;
xlabel('Time (mins)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Rotation (deg)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Rotation', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); 

% switch num_drums
%     case 8
%         ylim([55, 90]);
%     case 4
%         ylim([100, 135]);
%     case 2
%         ylim([110, 185]);
%     case 1
%         ylim([90, 190]);
%     otherwise
%         ylim([65, 85]);
% end

subplot(4,2,4);
plot(t(1:end-1), du(1:end-1), 'LineWidth', 2);
grid on;
xlabel('Time (mins)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Rate of change (deg/s)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Control drum speed', 'FontSize', font_size, 'FontWeight', font_weight);
%ylim([-1.5, 1.5]);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold

subplot(4,2,5);
plot(t(1:end-1), rho(1:end-1)*1e5, 'LineWidth', 2);
grid on;
xlabel('Time (mins)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Reactivity (PCM)', 'FontSize', font_size, 'FontWeight', font_weight);
title('Core Reactivity', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold


subplot(4,2,6);
plot(t(1:end-1), xHistory(10,1:end-1), 'LineWidth', 2); % Corrected to plot rows
grid on;
xlabel('Time (mins)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Temperature (K)', 'FontSize', font_size, 'FontWeight', font_weight);
ylim([870 910]);
title('Temperature of fuel', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); 


subplot(4,2,7);
plot(t(1:end-1), xHistory(11,1:end-1), 'LineWidth', 2);
grid on;
xlabel('Time (mins)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Temperature (K)', 'FontSize', font_size, 'FontWeight', font_weight);
ylim([ 870 900]);
title('Temperature of moderator', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold

subplot(4,2,8);
plot(t(1:end-1), xHistory(12,1:end-1), 'LineWidth', 2);
grid on;
xlabel('Time (mins)', 'FontSize', font_size, 'FontWeight', font_weight);
ylabel('Temperature (K)', 'FontSize', font_size, 'FontWeight', font_weight);
ylim([ 870 890]);
title('Temperature of coolant', 'FontSize', font_size, 'FontWeight', font_weight);
set(gca, 'FontSize', font_size, 'FontWeight', font_weight); % Set tick labels to bold

% Adjust the layout so labels don't overlap
set(gcf, 'Position', [5, 5, 2000, 4500]);

saveas(gcf, 'NMPC_power_simulation.png');


% Function to set parameters based on number of drums
function [Rho_d0, Reactivity_per_degree, u0] = setParameters(num_drums)
    switch num_drums
        case 8
            Rho_d0 = -0.033085599;
            Reactivity_per_degree = 26.11e-5;
            u0 = 77.56;
        case 4
            Rho_d0 = -0.033085599;
            Reactivity_per_degree = 16.11e-5;
            u0 = 108.5;
            %u0 = 125.56;
        case 2
            Rho_d0=-0.033085599+0.0074;
            %Rho_d0 = -0.033085599 + 0.0087;
            %Rho_d0 = -0.033085599 + 0.0078+0.0009;  %0.0006;
            Reactivity_per_degree = 7.33e-5;
            u0 = 165.5;
            %u0 = 174.84;
        case 1
            Rho_d0=-0.033085599+0.0071+0.0082;
            %Rho_d0 = -0.033085599 + 0.016; 
            %Rho_d0 = -0.033085599 + 0.0078+0.0082;
            Reactivity_per_degree = 2.77e-5;
            u0 = 170;
            %u0 = 174.84;
        otherwise
            Rho_d0 = -0.033085599;
            Reactivity_per_degree = 26.11e-5;
            u0 = 77.56;
    end
end

% Function to compute reactor dynamics
function dx = reactorDT0(x, u, Ts)
    
    % Define number of drums
    num_drums = 8;  % Define number of drums
    [Rho_d0, Reactivity_per_degree, ~] = setParameters(num_drums); 

    M = 5;  
    delta = Ts / M;
    xk1 = x;
    for ct = 1:M
        f1 = reactorCT0(xk1, u, Rho_d0, Reactivity_per_degree); 
        hxk1 = xk1 + delta * f1;
        f2 = reactorCT0(hxk1, u, Rho_d0, Reactivity_per_degree); 
        xk1 = xk1 + delta * (f1 + f2) / 2;
    end
    dx = xk1;
end

% Function for reactor continuous-time dynamics
function [dx, rho] = reactorCT0(x, u, Rho_d0, Reactivity_per_degree)
    % Constants and parameters
    Sig_x = 2.65e-22;
    yi = 0.061;
    yx = 0.002;
    lamda_x = 2.09e-5;
    lamda_I = 2.87e-5;
    Sum_f = 0.3358;

    % Coefficients and constants
    l = 1.68e-3;
    beta = 0.0048;
    beta_1 = 1.42481E-04; 
    beta_2 = 9.24281E-04;
    beta_3 = 7.79956E-04;
    beta_4 = 2.06583E-03;
    beta_5 = 6.71175E-04;
    beta_6 = 2.17806E-04;
    Lamda_1 = 1.272E-02;
    Lamda_2 = 3.174E-02;
    Lamda_3 = 1.160E-01;
    Lamda_4 = 3.110E-01;
    Lamda_5 = 1.400E+00;
    Lamda_6 = 3.870E+00;

    cp_f = 977;
    cp_m = 1697;
    cp_c = 5188.6;
    M_f = 2002;
    M_m = 11573;
    M_c = 500;
    mu_f = M_f * cp_f;
    mu_m = M_m * cp_m;
    mu_c = M_c * cp_c;
    f_f = 0.96;
    P_0 = 22e6; 
    Tf0 = 1105;
    Tm0 = 1087;
    T_in = 864;
    T_out = 1106;
    Tc0 = (T_in + T_out) / 2;
    K_fm = f_f * P_0 / (Tf0 - Tm0);
    K_mc = P_0 / (Tm0 - Tc0);
    M_dot = 1.75E+01;
    alpha_f = -2.875e-5;
    alpha_m = -3.696e-5;
    alpha_c = 0.0;
    X0 = 2.35496411413791e10;

    % State variables
    n_r = x(1); 
    Cr1 = x(2); 
    Cr2 = x(3); 
    Cr3 = x(4); 
    Cr4 = x(5);
    Cr5 = x(6);
    Cr6 = x(7);
    X = x(8); 
    I = x(9); 
    Tf = x(10); 
    Tm = x(11); 
    Tc = x(12); 

    % Reactivity and power calculations
    Rho_d1 = Rho_d0 + u * Reactivity_per_degree;
    %fprintf('Rho_d1: %.4f\n', Rho_d1);
    
    G = 3.2e-11;
    V = 400 * 200; 
    Pi = P_0 / (G * Sum_f * V);

    % Differential equations
    dx = zeros(12, 1);
    %rho = Rho_d1 + alpha_f * (Tf - Tf0) + alpha_c * (Tc - Tc0) + alpha_m * (Tm - Tm0) - Sig_x * (X - X0) / Sum_f;% 8,2,1 drums
    %rho     = Rho_d1 + alpha_f * (Tf - 900.42) + alpha_c * (Tc - 888.261) + alpha_m * (Tm - 898.261) - Sig_x * (X - X0) / Sum_f; %4drums

    dx(1) = (rho - beta) / l * n_r + beta_1 / l * Cr1 + beta_2 / l * Cr2 + beta_3 / l * Cr3 + beta_4 / l * Cr4 + beta_5 / l * Cr5 + beta_6 / l * Cr6;
    dx(2) = Lamda_1 * n_r - Lamda_1 * Cr1;
    dx(3) = Lamda_2 * n_r - Lamda_2 * Cr2;
    dx(4) = Lamda_3 * n_r - Lamda_3 * Cr3;
    dx(5) = Lamda_4 * n_r - Lamda_4 * Cr4;
    dx(6) = Lamda_5 * n_r - Lamda_5 * Cr5;
    dx(7) = Lamda_6 * n_r - Lamda_6 * Cr6;

    dx(8) = yx * Sum_f * Pi + lamda_I * I - Sig_x * X * Pi - lamda_x * X;
    dx(9) = yi * Sum_f * Pi - lamda_I * I;

    dx(10) = f_f * P_0 / mu_f * n_r - K_fm / mu_f * (Tf - Tc);
    dx(11) = (1 - f_f) * P_0 / mu_m * n_r + (K_fm * (Tf - Tm) - K_mc * (Tm - Tc)) / mu_m;
    dx(12) = K_mc * (Tm - Tc) / mu_c - 2 * M_dot * cp_c * (Tc - T_in) / mu_c;

    algebraic = [];
end

% Output function to be used by nmpc object 
function y = reactorOutputFcn(x, u, params, num_drums)
    y = x(1);
end

% Measurement function
function y = reactorMeasurementFcn(xk)
    y = xk(1);
end

% State function for extended Kalman filter
function xk1 = reactorStateFcn(xk, u)
    Ts = u(2); % Time step
    xk1 = reactorDT0(xk, u(1), Ts);
end
