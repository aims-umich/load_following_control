
This repository contains the design and implementation of a load-following control system for a microreactor. It includes three distinct control strategies:  

1. PID Control System 
2. Second-Order Super-Twisting Sliding Mode Control (SOSMC) 
3. Nonlinear Model Predictive Control (NMPC)

The microreactor model is based on the **nonlinear point kinetics equation**, coupled with:  

- A thermal-hydraulics model  
- A thermal neutron poison model 
- A reactivity model 

The reactor parameters are derived from the **Holos-Gen microreactor** developed by HolosGen LLC.  

The system is designed to track an **input power profile**, ensuring that the output power follows the desired trajectory with minimal error. The power profile can be dynamically adjusted based on demand.  

When executed with the appropriate **initial conditions** and **parameters**, the script generates the expected output, demonstrating the effectiveness of each control strategy.  
Each script is run separately to generate the power output. 
To run each script, ensure the number of drums correspond to the selected initial condition, and initial reactivity, which are indicated as commnents after each line. Only in the case of four drums that the initial reactivity differs but the initial reactivity works for the remaining. This is becasue the initial power is very low in the case of four drums. In addition, the the number of drums is selected in two different sections in the case of NMPC. The extra section is in the discretization function. 

Notes on running for MATLAB (this section must be removed after code optimization by Kamal):

For PID.m and STC.m:

- For 1, 2, and 8 drums: Change the following variables: `num_drums = 1/2/8`, use the correct value for the variable `pow` for 1/2/8 drum, and the correct value for the variable `x0` for 1/2/8 drum. Run the simulation.
- For 4 drums: Change the following variables: `num_drums = 4`, use the variable `pow` for 4 drums, and variable `x0` for 4 drums. Then go to line 343 in PID.m or line 367 in STC.m in the function `reactorDAE` and use the variable `rho` assigned for the 4 drum case and comment the other default one.

For NMPC.m
- For 1, 2, 4, and 8 drums: Change the following variables: `num_drums = 1/2/4/8`, use the variable `pow` for 1/2/4/8 drum, and variable `x0` for 1/2/4/8 drum. Run the simulation. **This simulation takes long time**.

Important Note: Please, note that Python scripts were not yet fully verified and optimized, so use with caution. 
