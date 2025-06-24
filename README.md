# Composite Hydrogen Tank Design for FCEVs 🚗💨

This repository presents the design, simulation, and analysis of a **Type III composite-wrapped hydrogen storage tank** for **Fuel Cell Electric Vehicles (FCEVs)**. The system is engineered to withstand a nominal pressure of **700 bar**, balancing **lightweight performance**, **structural safety**, and **manufacturing efficiency**.

> 🔧 Developed using MATLAB for structural simulation, this project integrates Classical Laminate Theory (CLT) and Tsai-Wu failure analysis for predicting ply-by-ply progressive failure in carbon fiber/epoxy laminates.

---

## 🔍 Project Overview

- **Goal**: Develop a modular, high-pressure hydrogen tank using lightweight composite materials to enhance safety, efficiency, and manufacturability.
- **Application**: Fuel Cell Electric Vehicles (FCEVs) – compatible with both sedans and SUVs.
- **Structure**:
  - **Liner**: Hydroformed 6061-T6 Aluminum
  - **Overwrap**: M60J Carbon Fiber / Aerospace Epoxy
- **Design Pressure**: 700 bar nominal

---

## 📊 Simulation & Methodology

- **Language**: MATLAB
- **Theories Used**:
  - Classical Laminate Theory (CLT)
  - Tsai-Wu Failure Criterion
  - Hygrothermal Stress Modeling

### 🔧 Key Features:
- Progressive **ply-by-ply failure simulation**
- **A-B-D stiffness matrix** generation for laminate response
- **Carpet plots** for strength ratio and deflection
- Residual stiffness modeling post-failure
- Integration of **thermal and moisture expansion** effects

---

## 📈 Results: Ply-by-Ply Stress-Strain Behavior

The figure below illustrates progressive ply failure in the composite laminate under increasing axial strain:

![image](https://github.com/user-attachments/assets/2772fcca-f9bf-4a55-b205-3c438822a11e)

---

### 🔎 Explanation:

- The **X-axis** shows the applied normal strain (εₓ), and the **Y-axis** shows the resulting axial stress (Nₓ).
- The simulation captured **progressive failure** of plies based on the **Tsai-Wu strength ratio**, and the laminate stiffness was reduced accordingly.
- Key events:
  - **First Ply Fail (90°)** at ~918 MPa: Indicates the most vulnerable orientation due to low transverse strength.
  - **±45° Ply Fail** at ~267 MPa: Fails under shear-dominated loading as axial strain increases.
  - **Last Ply Fail (0°)** at ~72 MPa: Final structural collapse once longitudinal plies exceed their limits.
- This behavior confirms **safe load-sharing** among orientations and highlights laminate failure sequence under multi-axial stress.

---

## 🛠️ Manufacturing Considerations

- **Hydroforming** used for the aluminum liner enables cost-effective mass production with excellent dimensional accuracy.
- **Filament winding** automates the composite overwrap for optimal fiber alignment and repeatability.

---

## ✅ Key Outcomes

- Verified safe composite layup under high-pressure conditions.
- Optimized ply orientation for strength, stiffness, and failure resistance.
- Demonstrated thermal and moisture stability through hygrothermal modeling.
- Achieved modular design concept for multi-cylinder integration in vehicles.

---

## 📁 Folder Structure

## 🤝 Let's Connect

If you're a recruiter or engineer interested in composites, energy storage, or structural simulation:

🔗 LinkedIn: (www.linkedin.com/in/baldeepsingh30)  
📫 Email: baldeep@lakeheadu.ca

📂 Full Project: [GitHub Repository Link](https://github.com/baldeepsingh30/composite-hydrogen-tank-design.git)


