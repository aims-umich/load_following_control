
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
To run each script, ensure the number of drums correspond to the selected initial condition, and initial reactivity, 



*For PID.m, STC.m, and NMPC.m*:

For 1, 2, 4, and 8 drums: Change the following variable: `num_drums = 1/2/4/8`
Run the simulation

*For NMPC.m:*
  
.**This simulation takes a long time**.

Important Note: Please, note that Python scripts were not yet fully verified and optimized, so use with caution. 
