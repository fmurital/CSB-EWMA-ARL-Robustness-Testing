# CSB-EWMA ARL₁ Robustness Testing

This repository contains code for testing the robustness of the CSB-EWMA control chart across different data-generating processes (normal, Laplace, uniform, and exponential).

## Features
- **Multiple Distributions**: Simulates and analyzes data from normal, Laplace, uniform, and exponential distributions.
- **ARL₁ Simulation**: Robust ARL₁ calculation under varying shifts (delta) for in-control ARL₀ targets (370 and 500).
- **Parallel Computing**: Utilizes parallel processing for large-scale simulations.
- **Dynamic UCL/LCL Limits**: Uses time-varying control limits for simulation.
  
## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/csb-ewma-arl1-robustness.git
