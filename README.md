# Vehicle-to-Grid (V2G) Bidirectional Power Flow Control

## Overview
This project simulates a **Vehicle-to-Grid (V2G)** system using MATLAB/Simulink. It demonstrates how Electric Vehicle (EV) batteries can function as energy storage units to stabilize grid frequency during peak load demand.

## Key Features
- **Bidirectional DC-DC Converter:** Controls power flow between the EV battery and the DC link.
- **Grid Synchronization:** Uses a Phase Locked Loop (PLL) to synchronize the inverter with the grid.
- **V2G & G2V Modes:**
  - **G2V (Grid-to-Vehicle):** Charges the battery when grid load is low.
  - **V2G (Vehicle-to-Grid):** Discharges battery power back to the grid during peak demand.

## Results
*(Upload your screenshot here! Click the image icon in the editor toolbar to add it)*
> The simulation shows the battery current reversing direction when the grid signals a peak load event, stabilizing the DC link voltage.

## Tech Stack
- **Software:** MATLAB & Simulink (R2023a)
- **Components:** IGBT Inverters, PID Controllers, LCL Filters.

## How to Run
1. Download the repository.
2. Open `V2G_Model.slx` in Simulink.
3. Run the simulation and view the scope results for Grid Voltage vs Battery Current.
