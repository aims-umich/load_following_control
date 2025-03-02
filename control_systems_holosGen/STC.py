#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: abkamal
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d

# Supertwisting Controller Class
class STCController:
    def __init__(self, lambda_, W, T_C, max_value, min_value, max_rate, u0):
        self.lambda_ = lambda_
        self.W = W
        self.T_C = T_C
        self.max_value = max_value
        self.min_value = min_value
        self.max_rate = max_rate
        
        self.w = u0
        self.err_prev = 0
        self.deriv_prev = 0
        self.t_prev = 0
        self.command_sat_prev = u0
        self.command_prev = u0

    def update(self, t, measurement, setpoint):
        err = setpoint - measurement
        T = t - self.t_prev
        self.t_prev = t
        
        deriv_filt = (err - self.err_prev + self.T_C * self.deriv_prev) / (T + self.T_C)
        self.err_prev = err
        self.deriv_prev = deriv_filt
        
        s = err + deriv_filt
        self.w += self.W * np.sign(s) * T
        command = self.lambda_ * np.sqrt(np.abs(s)) * np.sign(s) + self.w
        self.command_prev = command
        
        command_sat = np.clip(command, self.min_value, self.max_value)
        if command_sat > self.command_sat_prev + self.max_rate * T:
            command_sat = self.command_sat_prev + self.max_rate * T
        elif command_sat < self.command_sat_prev - self.max_rate * T:
            command_sat = self.command_sat_prev - self.max_rate * T

        self.command_sat_prev = command_sat

        return command_sat

# Simulation parameters
dt = 0.1
T = 2000
n0 = 1
time = np.arange(0, T + dt, dt)
nt = len(time)
ref = np.zeros(nt)

num_drums = 1

if num_drums == 8:
    Rho_d0 = -0.033085599
    Reactivity_per_degree = 26.11e-5
    u0 = 77.56
    W = 0.06
elif num_drums == 4:
    Rho_d0 = -0.033085599
    Reactivity_per_degree = 16.11e-5
    u0 = 125.56
    W = 0.1
elif num_drums == 2:
    Rho_d0 = -0.033085599 + 0.0073
    Reactivity_per_degree = 7.33e-5
    u0 = 174.84
    W = 0.2
elif num_drums == 1:
    Rho_d0 = -0.033085599 + 0.0071 + 0.0082
    Reactivity_per_degree = 2.77e-5
    u0 = 174.84
    W = 0.49
else:
    Rho_d0 = -0.033085599
    Reactivity_per_degree = 26.11e-5
    u0 = 77.56
    W = 0.04

time_point = np.array([0, 20, 30, 50, 60, 80, 90, 110, 130, 200]) * 10
#pow = np.array([1, 1, 0.4, 0.4, 1, 1, 0.4, 0.4, 1, 1])

#time_point=np.array([0, 20, 30, 50, 60, 80, 90, 110, 130, 150])
#time_point=np.array([0, 10, 40, 70, 80, 100, 110, 130, 140, 150])
#time_point=np.array([0, 30, 40, 50, 60, 70, 90, 100, 140, 150])
#time_point=np.array([0, 20, 30, 50, 60, 80, 90, 110, 130, 200])*10;
#pow=np.array([0.8, 0.8, 0.4, 0.4, 0.8, 0.8, 0.4, 0.4, 0.8, 0.8])
pow=np.array([1, 1, 0.4, 0.4, 1, 1, 0.4, 0.4, 1, 1])
#pow=np.array([0.3, 1, 0.4, 0.4, 0.6, 0.6, 0.8, 0.8, 1, 1])
#pow=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
#pow=np.array([0.9, 0.9, 0.3, 0.3, 0.9, 0.9, 0.3, 0.3, 0.9, 0.9])
#pow=np.array([0.7, 0.7, 0.4, 0.4, 0.4, 0.4, 0.6, 0.6, 0.7, 0.7])
#pow=np.array([0.5, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

ref_old = pow[0]
for it in range(nt):
    if it > 0:
        time[it] = time[it - 1] + dt
        ref_old = ref[it - 1]
    ref[it] = ref_old
    for ii in range(len(time_point) - 1):
        if time_point[ii] <= time[it] <= time_point[ii + 1]:
            frac1 = (time_point[ii + 1] - time[it]) / (time_point[ii + 1] - time_point[ii])
            frac2 = 1.0 - frac1
            ref[it] = frac1 * pow[ii] + frac2 * pow[ii + 1]
            break

yref = ref * n0

# STC tuning parameters
lambda_ = 0.001
T_C = 0.5
max_value = 180
min_value = 0
max_rate = 0.5

# Initial conditions
Sig_x = 2.65e-22
yi = 0.061
yx = 0.002
lamda_x = 2.09e-5
lamda_I = 2.87e-5
Sum_f = 0.3358
G = 3.2e-11
V = 400 * 200
P_0 = 22e6
Pi = P_0 / (G * Sum_f * V)

Xe0 = (yi + yx) * Sum_f * Pi / (lamda_x + Sig_x * Pi)
I0 = yi * Sum_f * Pi / lamda_I

x0 = np.array([pow[0], pow[0], pow[0], pow[0], pow[0], pow[0], pow[0], I0, Xe0, 900.42, 898.28, 888.261])

# Instantiate the STC controller
stc_controller = STCController(lambda_, W, T_C, max_value, min_value, max_rate, u0)

# Reactor dynamics and STC controller integration
u_stc = np.zeros(nt)
x_stc = np.zeros((nt, len(x0)))
rho = np.zeros(nt)  # Initialize rho as an array
x_stc[0] = x0
u_stc[0] = u0

def reactorDAE(t, x, u, Rho_d0, Reactivity_per_degree, Xe0, I0, Pi):
    # Reactor parameters
    Sig_x = 2.65e-22
    yi = 0.061
    yx = 0.002
    lamda_x = 2.09e-5
    lamda_I = 2.87e-5
    Sum_f = 0.3358
    l = 1.68e-3
    beta = 0.0048
    beta_1 = 1.42481e-04
    beta_2 = 9.24281e-04
    beta_3 = 7.79956e-04
    beta_4 = 2.06583e-03
    beta_5 = 6.71175e-04
    beta_6 = 2.17806e-04
    Lamda_1 = 1.272e-02
    Lamda_2 = 3.174e-02
    Lamda_3 = 1.160e-01
    Lamda_4 = 3.110e-01
    Lamda_5 = 1.400e+00
    Lamda_6 = 3.870e+00
    cp_f = 977
    cp_m = 1697
    cp_c = 5188.6
    M_f = 2002
    M_m = 11573
    M_c = 500
    mu_f = M_f * cp_f
    mu_m = M_m * cp_m
    mu_c = M_c * cp_c
    f_f = 0.96
    P_0 = 22e6
    Tf0 = 1105
    Tm0 = 1087
    T_in = 864
    T_out = 1106
    Tc0 = (T_in + T_out) / 2
    K_fm = f_f * P_0 / (Tf0 - Tm0)
    K_mc = P_0 / (Tm0 - Tc0)
    M_dot = 1.75e01
    alpha_f = -2.875e-5
    alpha_m = -3.696e-5
    alpha_c = 0.0
    
    # State variables
    n_r, Cr1, Cr2, Cr3, Cr4, Cr5, Cr6, X, I, Tf, Tm, Tc = x
    Rho_d1 = Rho_d0 + u * Reactivity_per_degree
    
    # ODEs
    dx = np.zeros(12)
    rho = Rho_d1 + alpha_f * (Tf - Tf0) + alpha_c * (Tc - Tc0) + alpha_m * (Tm - Tm0) - Sig_x * (X - Xe0) / Sum_f
    
    # Kinetics equations with six-delayed neutron groups
    dx[0] = (rho - beta) / l * n_r + sum([
        beta_1 / l * Cr1, 
        beta_2 / l * Cr2, 
        beta_3 / l * Cr3, 
        beta_4 / l * Cr4, 
        beta_5 / l * Cr5, 
        beta_6 / l * Cr6
    ])
    dx[1] = Lamda_1 * n_r - Lamda_1 * Cr1
    dx[2] = Lamda_2 * n_r - Lamda_2 * Cr2
    dx[3] = Lamda_3 * n_r - Lamda_3 * Cr3
    dx[4] = Lamda_4 * n_r - Lamda_4 * Cr4
    dx[5] = Lamda_5 * n_r - Lamda_5 * Cr5
    dx[6] = Lamda_6 * n_r - Lamda_6 * Cr6
    
    # Xenon and Iodine dynamics
    dx[7] = yx * Sum_f * Pi + lamda_I * I - Sig_x * X * Pi - lamda_x * X
    dx[8] = yi * Sum_f * Pi - lamda_I * I
    
    # Thermal–hydraulics model of the reactor core
    dx[9] = f_f * P_0 / mu_f * n_r - K_fm / mu_f * (Tf - Tc)
    dx[10] = (1 - f_f) * P_0 / mu_m * n_r + (K_fm * (Tf - Tm) - K_mc * (Tm - Tc)) / mu_m
    dx[11] = K_mc * (Tm - Tc) / mu_c - 2 * M_dot * cp_c * (Tc - T_in) / mu_c
    
    return dx, rho

for i in range(nt - 1):
    dx, rho[i] = reactorDAE(0, x_stc[i], u_stc[i], Rho_d0, Reactivity_per_degree, Xe0, I0, Pi)
    x_stc[i + 1] = x_stc[i] + dx * dt
    u_stc[i + 1] = stc_controller.update(time[i], x_stc[i + 1, 0], ref[i + 1])

# Final reactivity calculation
dx, rho[nt-1] = reactorDAE(0, x_stc[nt-1], u_stc[nt-1], Rho_d0, Reactivity_per_degree, Xe0, I0, Pi)

t = time
du_stc = np.concatenate(([0], np.diff(u_stc) / dt))

# Low-pass filter implementation (simple moving average)
window_size = 30
du_stc_smooth = np.copy(du_stc)
for i in range(window_size, len(du_stc) - window_size):
    du_stc_smooth[i] = np.mean(du_stc[i - window_size:i + window_size])

# Interpolation
interp_ref = interp1d(time, ref, kind='linear', fill_value="extrapolate")
r = interp_ref(t)
y = x_stc[:, 0]

# Error analysis
error = y - r
downsample_factor = 10
t_ds = t[::downsample_factor]
r_ds = r[::downsample_factor]
y_ds = y[::downsample_factor]

MAE = np.mean(np.abs(r_ds - y_ds))
control_effort = np.linalg.norm(du_stc)

print(f'Mean Absolute Error (MAE): {MAE:.8f}')
print(f'Control Effort: {control_effort:.2f}')

# Plot results
plt.figure(figsize=(10, 15))
gs = gridspec.GridSpec(4, 2, height_ratios=[1, 1, 1, 1])

# First subplot: occupies the top row (first and second positions)
ax1 = plt.subplot(gs[0, :])  # Span the entire top row
ax1.plot(t, y * 100, linewidth=2)
ax1.plot(t, r * 100, '--', linewidth=2)
ax1.grid(True)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Power (%)')
ax1.set_title('STC microreactor control system core power simulation with 1 control drum')
ax1.set_ylim([0, 200])
ax1.legend(['Actual power', 'Desired power'])

# Second subplot: Rotation
ax2 = plt.subplot(gs[1, 0])
ax2.plot(time, u_stc, linewidth=2)
ax2.grid(True)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Rotation (deg)')
ax2.set_title('Rotation')
if num_drums == 8:
    ax2.set_ylim([65, 85])
elif num_drums == 4:
    ax2.set_ylim([110, 135])
elif num_drums == 2:
    ax2.set_ylim([120, 185])
elif num_drums == 1:
    ax2.set_ylim([120, 190])
else:
    ax2.set_ylim([65, 85])

# Third subplot: Rate of change
ax3 = plt.subplot(gs[1, 1])
ax3.plot(time, du_stc_smooth, linewidth=2)
ax3.grid(True)
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Rate of change (deg/s)')
ax3.set_title('Rate of change')
ax3.set_ylim([-1, 1])

# Fourth subplot: Reactivity over Time
ax4 = plt.subplot(gs[2, 0])
ax4.plot(t[:-1], rho[:-1] * 1e5, linewidth=2)
ax4.grid(True)
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('Reactivity (PCM)')
ax4.set_title('Reactivity over Time')

# Fifth subplot: Temperature of fuel
ax5 = plt.subplot(gs[2, 1])
ax5.plot(t, x_stc[:, 9], linewidth=2)
ax5.grid(True)
ax5.set_xlabel('Time (s)')
ax5.set_ylabel('Temperature (K)')
ax5.set_title('Temperature of fuel')

# Sixth subplot: Temperature of moderator
ax6 = plt.subplot(gs[3, 0])
ax6.plot(t, x_stc[:, 10], linewidth=2)
ax6.grid(True)
ax6.set_xlabel('Time (s)')
ax6.set_ylabel('Temperature(K)')
ax6.set_title('Temperature of moderator')

# Seventh subplot: Temperature of coolant
ax7 = plt.subplot(gs[3, 1])
ax7.plot(t, x_stc[:, 11], linewidth=2)
ax7.grid(True)
ax7.set_xlabel('Time (s)')
ax7.set_ylabel('Temperature (K)')
ax7.set_title('Temperature of coolant')

plt.tight_layout()
plt.savefig('STC_power_simulation.png')
plt.show()

