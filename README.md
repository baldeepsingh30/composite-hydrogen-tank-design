# composite-hydrogen-tank-design
This repository contains the complete design, simulation, and analysis of a Type III composite-wrapped hydrogen storage tank for fuel cell electric vehicles (FCEVs), including MATLAB scripts for laminate analysis and failure modeling using Classical Laminate Theory and Tsai-Wu criteria.

## 🔍 Project Overview

- **Goal**: Develop a modular, high-pressure hydrogen tank using lightweight composite materials to enhance safety, efficiency, and manufacturability.
- **Application**: Fuel Cell Electric Vehicles (FCEVs) – suitable for sedans and SUVs.
- **Structure**: 
  - Liner: 6061-T6 aluminum, hydroformed
  - Overwrap: Carbon Fiber/Epoxy (M60J)
- **Pressure Rating**: 700 bar nominal

---

## 🧪 Simulation Details

- **Simulation Language**: MATLAB
- **Theory Applied**:
  - Classical Laminate Theory (CLT)
  - Tsai-Wu Failure Criterion
  - Hygrothermal stress analysis
- **Key Features**:
  - Ply-by-ply progressive failure simulation
  - A-B-D stiffness matrix generation
  - Carpet plots for strength and deflection
  - Residual stiffness modeling after ply failure

---
### ✅ Additional Insight:

This ply-by-ply failure simulation of the carbon fiber/epoxy laminate confirms that the composite shell can **withstand internal pressures exceeding 700 bar**, as the **first ply failure occurs at approximately 918 MPa** of axial stress—well above the design threshold.

To account for the **load-sharing with the aluminum liner**, a compensating equation was developed:

![Load Sharing Equation](https://quicklatex.com/cache3/c2/ql_b70fa99730189ddb90d2fb2e926413c2_l3.png)

Where:
- `σ_total`: Total internal pressure (converted to equivalent axial stress)
- `σ_aluminum`: Hoop stress taken by the aluminum liner
- `t_al`, `t_composite`, `t_total`: Thicknesses of aluminum, composite, and total wall

This formulation allowed for accurate estimation of the pressure **absorbed solely by the composite laminate**, validating its performance under realistic structural conditions.

---

## 🛠️ Manufacturing Considerations

- **Liner**: Hydroforming technique for cost-effective shaping
- **Composite Overwrap**: Automated filament winding for performance and scalability

---

## 📈 Outcomes

- Verified composite layup under operational loads
- Optimized fiber orientation for strength and weight
- Designed a compact, modular layout for multi-cylinder integration
- Validated results through stress-strain tracking and thermal simulations

---

## 📁 Folder Structure
