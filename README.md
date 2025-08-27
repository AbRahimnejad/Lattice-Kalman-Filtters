# Lattice Kalman Filter – MATLAB Implementation

## 📌 Overview
This package demonstrates the **Lattice Kalman Filter (LKF)** applied to a nonlinear **Electro-Hydrostatic Actuator (EHA)** system.  
The filter leverages **Gaussianized Korobov lattice points** with a **Cranley–Patterson random shift** to approximate integrals in the prediction and update steps.  

The workflow:
1. Simulate the nonlinear plant with process and measurement noise.  
2. Run the LKF to estimate the states.  
3. Report RMSE and generate comparison plots.

---

## 📂 File Structure
- **`main.m`** – Entry point script (simulation + plotting).  
- **`lkf.m`** – Implements the Lattice Kalman Filter.  
- **`mod1shift.m`** – Utility for Cranley–Patterson random shifting of lattice points.

---

## ⚙️ Requirements
- MATLAB R2020a or later (earlier versions work if `erfinv` is available).  
- No additional toolboxes required.

---

## 🚀 How to Run
1. Place all three files in the same working folder:  
   - `main.m`  
   - `lkf.m`  
   - `mod1shift.m`  
2. Open MATLAB and set this folder as the **Current Folder**.  
3. Run:
   ```matlab
   main

## 📖 Reference

For the theoretical foundation and original introduction of the Lattice Kalman Filter, please refer to:

**A. Rahimnejad, S. A. Gadsden, and M. Al-Shabi,**  
*“Lattice Kalman Filters,”* IEEE Signal Processing Letters, vol. 28, pp. 1355–1359, 2021.  
[https://doi.org/10.1109/LSP.2021.3089935](https://doi.org/10.1109/LSP.2021.3089935)

