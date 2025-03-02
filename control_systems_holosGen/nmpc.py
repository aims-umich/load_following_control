import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import time

# Define constants and initial parameters
nx = 12
ny = 1
nu = 1
dt = 0.5  # Increased time step for faster simulation
Ts = 0.05  # Adjusted sampling time

# MPC parameters
PredictionHorizon = 5  # Reduced horizon for faster computation
ControlHorizon = 2  # Reduced control horizon

# Define number of drums
num_drums = 8  # Define number of drums

# Simulation parameters
T = 6000
sim_time = np.arange(0, T + dt, dt)  # Renamed from `time` to `sim_time`
nt = len(sim_time)
ref = np.zeros(nt)

# Function to set parameters based on number of drums
def setParameters(num_drums):
    if num_drums == 8:
        Rho_d0 = -0.033085599
        Reactivity_per_degree = 26.11e-5
        u0 = 77.56
    elif num_drums == 4:
        Rho_d0 = -0.033085599
        Reactivity_per_degree = 16.11e-5
        u0 = 125.56
    elif num_drums == 2:
        Rho_d0 = -0.033085599 + 0.0073
        Reactivity_per_degree = 7.33e-5
        u0 = 177.84
    elif num_drums == 1:
        Rho_d0 = -0.033085599 + 0.0078 + 0.0082
        Reactivity_per_degree = 2.77e-5
        u0 = 177.84
    else:
        Rho_d0 = -0.033085599
        Reactivity_per_degree = 26.11e-5
        u0 = 77.56
    return Rho_d0, Reactivity_per_degree, u0

Rho_d0, Reactivity_per_degree, u0 = setParameters(num_drums)

# Time-varying reference for simulation
time_point = np.array([0, 20, 30, 50, 60, 80, 90, 110, 130, 200]) * 30
pow_point = np.array([1, 1, 0.4, 0.4, 1, 1, 0.4, 0.4, 1, 1])
#pow_point=np.array([0.7, 0.7, 1, 1, 0.9, 0.9, 1, 1, 0.5, 0.5]) 
#pow_point=np.array([0.7, 0.7, 0.4, 0.4, 0.4, 0.8, 0.8, 0.8, 1, 1]) 
ref_old = pow_point[0]

# Generate reference trajectory
for it in range(nt):
    if it > 0:
        ref_old = ref[it - 1]
    ref[it] = ref_old
    for ii in range(len(time_point) - 1):
        if time_point[ii] <= sim_time[it] <= time_point[ii + 1]:
            frac1 = (time_point[ii + 1] - sim_time[it]) / (time_point[ii + 1] - time_point[ii])
            frac2 = 1.0 - frac1
            ref[it] = frac1 * pow_point[ii] + frac2 * pow_point[ii + 1]
            break

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



# Initialize state and EKF

x0 = np.array([pow_point[0], pow_point[0], pow_point[0], pow_point[0], pow_point[0], pow_point[0], pow_point[0], Xe0, I0, 900.42, 898.28, 888.261])
x = x0
y = x0[0]
mv = u0

# EKF initializations
P = np.eye(nx)  # Initial covariance matrix
Q = 0.01 * np.eye(nx)  # Process noise covariance
R = 0.1  # Measurement noise covariance
H = np.array([[1] + [0] * (nx - 1)])  # Measurement matrix and it ensures the first element is 1 and the ramiining are 0

# Initialize simulation variables
Duration = 6000
xHistory = np.zeros((nx, len(sim_time)))
xkHistory = np.zeros((nx, len(sim_time)))
MV = np.zeros(len(sim_time))

xHistory[:, 0] = x0
xkHistory[:, 0] = x0
MV[0] = u0

"""
The function numericalJacobian computes the numerical approximation
of the Jacobian matrix for a given function fun
with respect to the state vector x. This is done using finite difference approximation.
"""

# Numerical Jacobian calculation
def numericalJacobian(fun, x, u, Ts):
    n = len(x) # # Get the size of the state vector
    epsilon = 1e-5 # # Small perturbation for numerical differentiation/finite difference approximation
    J = np.zeros((n, n)) # Initialize the Jacobian matrix as an nxn zero matrix
    fx = fun(x, np.array([u, Ts]))
    for i in range(n): # Loop through each state variable
        x_eps = x.copy() # Copy the state vector
        x_eps[i] += epsilon # Perturb the i-th state slightly
        fx_eps = fun(x_eps, np.array([u, Ts]))  # Evaluate function at perturbed state
        J[:, i] = (fx_eps - fx) / epsilon # Compute numerical derivative
    return J

# Function to compute reactor dynamics
def reactorCT0(x, u, Rho_d0, Reactivity_per_degree):
    # Constants and parameters
    Sig_x = 2.65e-22
    yi = 0.061
    yx = 0.002
    lamda_x = 2.09e-5
    lamda_I = 2.87e-5
    Sum_f = 0.3358
    l = 1.68e-3
    beta = 0.0048
    beta_1 = 1.42481e-4
    beta_2 = 9.24281e-4
    beta_3 = 7.79956e-4
    beta_4 = 2.06583e-3
    beta_5 = 6.71175e-4
    beta_6 = 2.17806e-4
    Lamda_1 = 1.272e-2
    Lamda_2 = 3.174e-2
    Lamda_3 = 1.16e-1
    Lamda_4 = 3.11e-1
    Lamda_5 = 1.4e0
    Lamda_6 = 3.87e0

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
    M_dot = 1.75e1
    alpha_f = -2.875e-5
    alpha_m = -3.696e-5
    alpha_c = 0.0
    X0 = 2.35496411413791e10

    # State variables
    n_r, Cr1, Cr2, Cr3, Cr4, Cr5, Cr6, X, I, Tf, Tm, Tc = x

    # Reactivity and power calculations
    Rho_d1 = Rho_d0 + u * Reactivity_per_degree

    G = 3.2e-11
    V = 400 * 200
    Pi = P_0 / (G * Sum_f * V)

    # Differential equations
    dx = np.zeros(12)
    rho = (Rho_d1 + alpha_f * (Tf - Tf0) + alpha_c * (Tc - Tc0) +
           alpha_m * (Tm - Tm0) - Sig_x * (X - X0) / Sum_f)

    dx[0] = (rho - beta) / l * n_r + beta_1 / l * Cr1 + beta_2 / l * Cr2 + beta_3 / l * Cr3 + beta_4 / l * Cr4 + beta_5 / l * Cr5 + beta_6 / l * Cr6
    dx[1] = Lamda_1 * n_r - Lamda_1 * Cr1
    dx[2] = Lamda_2 * n_r - Lamda_2 * Cr2
    dx[3] = Lamda_3 * n_r - Lamda_3 * Cr3
    dx[4] = Lamda_4 * n_r - Lamda_4 * Cr4
    dx[5] = Lamda_5 * n_r - Lamda_5 * Cr5
    dx[6] = Lamda_6 * n_r - Lamda_6 * Cr6

    dx[7] = yx * Sum_f * Pi + lamda_I * I - Sig_x * X * Pi - lamda_x * X
    dx[8] = yi * Sum_f * Pi - lamda_I * I

    dx[9] = f_f * P_0 / mu_f * n_r - K_fm / mu_f * (Tf - Tc)
    dx[10] = ((1 - f_f) * P_0 / mu_m * n_r +
              (K_fm * (Tf - Tm) - K_mc * (Tm - Tc)) / mu_m)
    dx[11] = (K_mc * (Tm - Tc) / mu_c - 2 * M_dot * cp_c * (Tc - T_in) / mu_c)

    return dx

# Function to compute discretized reactor dynamics using Heun's method  
def reactorDT0(x, u, Ts):
    num_drums = 8  # Define number of drums
    Rho_d0, Reactivity_per_degree, _ = setParameters(num_drums)  # Use default for dx calculation

    M = 5
    delta = Ts / M
    xk1 = x
    for _ in range(M):
        f1 = reactorCT0(xk1, u, Rho_d0, Reactivity_per_degree)
        hxk1 = xk1 + delta * f1
        f2 = reactorCT0(hxk1, u, Rho_d0, Reactivity_per_degree)
        xk1 = xk1 + delta * (f1 + f2) / 2
    return xk1

# Output function to be used by NMPC object 
def reactorOutputFcn(x, u, params, num_drums):
    return x[0]

# State function for extended Kalman filter
def reactorStateFcn(xk, u):
    mv = u[0]  # Manipulated variable
    Ts = u[1]  # Time step
    return reactorDT0(xk, mv, Ts)

# Manual NMPC implementation
def manualNMPC(xk, mv, ref, PredictionHorizon, ControlHorizon, Ts, nx, nu, dt):
    Q = 1e7  # Weight on output variable
    R = 0.1  # Weight on manipulated variable
    Rd = 1  # Weight on rate of change of manipulated variable

    u_opt = np.ones(ControlHorizon) * mv  # Initial guess for control inputs
    u_lb = 0  # Lower bound on manipulated variable
    u_ub = 180  # Upper bound on manipulated variable
    du_lb = -1 * dt  # Lower bound on rate of change of manipulated variable
    du_ub = 1 * dt  # Upper bound on rate of change of manipulated variable

    # Objective function
    def objFun(u):
        return nmpcObjective(u, xk, ref, PredictionHorizon, ControlHorizon, Q, R, Rd, Ts, nx, nu)

    # Constraints
    constraints = nmpcConstraints(u_opt, mv, ControlHorizon, u_lb, u_ub, du_lb, du_ub)
    
# uses SLSQP- sequential least squares quadratic programming. 
# bounds=[(u_lb, u_ub)] * ControlHorizon ---to ensure 

    """
    Ensures Each Control Input is Constrained
    Since we are optimizing a sequence of control inputs over ControlHorizon, each one needs bounds.
    Prevents scipy.optimize.minimize from Applying Bounds to Only the First Input
    Without this, the optimizer may only restrict the first input, allowing later inputs to be unbounded.
    Works for Any Control Horizon Length
    If ControlHorizon = 10, we automatically apply bounds to all 10 control inputs
    """ 

    result = minimize(objFun, u_opt, method='SLSQP', bounds=[(u_lb, u_ub)] * ControlHorizon, constraints=constraints)
    u_opt = result.x
    return u_opt[0] # obtain the first control input in the sequence and apply it to the system. 

# NMPC objective function
def nmpcObjective(u, xk, ref, PredictionHorizon, ControlHorizon, Q, R, Rd, Ts, nx, nu):
    J = 0
    x = xk.copy()

    for k in range(PredictionHorizon):
        x = reactorDT0(x, u[min(k, ControlHorizon-1)], Ts) #min(k, ControlHorizon-1) ensures that: If k < ControlHorizon, 
        #it uses the actual control input u[k]., If k ≥ ControlHorizon, it uses the last available control input u[ControlHorizon-1].
        y = reactorOutputFcn(x, [], [], 8)
        J += Q * (y - ref)**2
        if k < ControlHorizon:
            J += R * u[k]**2
            if k > 0:
                J += Rd * (u[k] - u[k-1])**2
    return J

# NMPC constraints
def nmpcConstraints(u, mv, ControlHorizon, u_lb, u_ub, du_lb, du_ub):
    cons = []

    for k in range(ControlHorizon):
        cons.append({'type': 'ineq', 'fun': lambda u, k=k: u[k] - u_lb})
        cons.append({'type': 'ineq', 'fun': lambda u, k=k: u_ub - u[k]})

    for k in range(1, ControlHorizon):
        cons.append({'type': 'ineq', 'fun': lambda u, k=k: u[k] - u[k-1] - du_lb})
        cons.append({'type': 'ineq', 'fun': lambda u, k=k: du_ub - (u[k] - u[k-1])})

    cons.append({'type': 'ineq', 'fun': lambda u: u[0] - mv - du_lb})
    cons.append({'type': 'ineq', 'fun': lambda u: du_ub - (u[0] - mv)})

    return cons

# Main simulation loop
start_time = time.time()
for ct in range(1, len(sim_time)):
    # EKF prediction step
    
    """
    F is the Jacobian matrix of the state transition function.
    P is the current covariance matrix.
    F.T is the transpose of F.
    Q is the process noise covariance matrix.
    
    the mathematical expression for the Kalman filter 
    Ppred= F*P*F' + Q
    K=Ppred* H'*inv(Hpred*H'+R)
    Xk=Xpred+K(y-H*xpred)
    P=(1-KH)Ppred
    
    
    """
    xk_pred = reactorDT0(x, mv, Ts)
    F = numericalJacobian(reactorStateFcn, x, mv, Ts)
    P_pred = F @ P @ F.T + Q
    
    # EKF update step
    K = P_pred @ H.T @ np.linalg.inv(H @ P_pred @ H.T + R)
    xk = xk_pred + K @ (y - H @ xk_pred)
    P = (np.eye(nx) - K @ H) @ P_pred
    
    # Manual NMPC control step
    mv = manualNMPC(xk, mv, ref[ct], PredictionHorizon, ControlHorizon, Ts, nx, nu, dt)

    # Update state
    x = reactorDT0(x, mv, Ts)
    y = x[0]

    # Store history
    xHistory[:, ct] = x
    xkHistory[:, ct] = xk
    MV[ct] = mv

# End time for simulation
end_time = time.time()
print(f"Simulation time: {end_time - start_time} seconds")

# Post-processing: compute and display MAE
error = xHistory[0, :] - ref
MAE = np.mean(np.abs(error))
print(f'Mean Absolute Error (MAE): {MAE:.8f}')

# Apply a moving average filter to smooth the rate of change of MV manually
windowSize = 10
du = np.diff(MV) / dt
du_smooth = np.convolve(du, np.ones(windowSize)/windowSize, mode='same')

# Plotting results
t = np.arange(0, Duration + dt, dt)
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(t, xHistory[0, :] * 100, label='Actual power', linewidth=2)
plt.plot(t, ref * 100, '--', label='Desired power', linewidth=2)
plt.grid(True)
plt.xlabel('Time (s)')
plt.ylabel('Power (%)')
plt.title('NMPC Microreactor System Core Power Simulation with 1 control drum')
plt.ylim([0, 200])
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(t, MV, linewidth=2)
plt.ylabel('Manipulated Variable')
plt.xlabel('Time (s)')
#plt.ylim([110, 185])

plt.subplot(3, 1, 3)
plt.plot(t[1:], du_smooth, linewidth=2)
plt.ylim([-1.5, 1.5])
plt.ylabel('Rate of Change of MV')
plt.xlabel('Time (s)')

plt.tight_layout()
plt.show()
