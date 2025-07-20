# Gravitational-Waveform-and-Orbit
Simulation of Gravitational Waves (Waveform and Orbits) from a Compact Binary Coalescence (CBC) in Newtonian Approximation




This repository contains code to simulate the inspiral of a compact binary system under gravitational wave emission, using the quadrupole approximation. The goal is to visualize both the gravitational waveform and the orbital evolution as the system loses energy and angular momentum.

## Features

- **Waveform and Orbit Plots (`waveform_and_orbit.ipynb`)**  
  Jupyter notebook that computes and plots:
  - The gravitational wave strain $h_+(t)$ for a compact binary on an eccentric orbit.
  - The evolution of the orbital trajectory as $e(t)$ (eccentricity) and  $p(t)$ (semi-latus rectum) decrease over time due to gravitational wave emission.

- **Combined Animation (`BNS_simulation.py`)**  
  Python script that animates both the gravitational waveform and the orbital motion in a synchronized view.  
  It shows how the waveform and the orbit change together, illustrating the circularization of the orbit and the decrease of the orbital radius as the system inspirals.

## Requirements

- Python 3.x
- `numpy`
- `matplotlib`
- `scipy`
- `jupyter` (for running the notebook)

Install dependencies with:
```bash
pip install numpy matplotlib scipy jupyter
````

## How to Run

1. **Visualize Waveform and Orbit in Jupyter Notebook:**

   ```bash
   jupyter notebook waveform_and_orbit.ipynb
   ```

2. **Run the Combined Animation:**

   ```bash
   python BNS_simulation.py
   ```

## Physical Context

The orbital evolution is computed using the leading-order post-Newtonian quadrupole formulas, which describe how gravitational radiation reaction causes the orbit to shrink and circularize over time.
This approach allows us to visualize:

* The orbital decay and circularization.
* The corresponding evolution of the gravitational waveform.

## License

MIT License — free to use, modify, and share.



