# Kirkwood
## Kirkwood gap simulation
## Introduction
The distribution of the semi-major axes of the orbits of asteroids in the inner solar system isn't homegeneous. Indeed, there are gaps and peaks, unevenly distributed. 
D. Kirkwood first noticed them in 1866 and suggested the gaps were due to orbital resonances with Jupiter.

## Simulation
### Hypothesis
The aim of the simulation is to reproduce the Kirkwood gaps from an initial distribution of asteroids from the main belt [1.8 to 4.2 AU, roughly], under the influence of the Sun, Jupiter and a few other major bodies (Mars, Saturn).
Most asteroids have a small inclination which is neglected here to reduce computational time.
Collisions between asteroids are rare events, as their average distance is ~ 10<sup>6</sup> km. Thus, the gravitational attraction from the main bodies will be the only effect accounted.

### Code
#### Parameters
Main parameters are set in `kirkwood.h` :
- `ASTEROID_NUMBER`   number of asteroids generated (won't work if it exceeds ~15000)
- `JUPITER_REV`       length of the simulation (in number of Jupiter's revolution around the Sun)
- `STEP_TIME`         works well enough up to ~0.2 year

#### Integrator
The numerical scheme used for the simulation is a 4th-order Yoshida integrator, which is derived from a leapfrog integrator.
It ensures the total energy is conserved through time.

#### Input files
`asteroids.dat` is usually generated at the beginning of the simulation if `RESET_ASTEROID_INPUT` is True.
The main bodies standard gravitational parameter and initial state vectors are stored within the other files. One can switch between:
- `jupiter.dat` (Sun & Jupiter)
- `mars.dat`    (same + Mars)
- `saturn.dat`  (same + Saturn)

#### Output files
Each simulation provides output data files:
- state vectors for each body (main bodies included) at each `STEP_WRITING` time, if `SAVE_BOOL` is True (although it increases significantly the simulation time)
- orbital elements (semi-major axis, eccentricity and true longitude of periapsis) at the initial and final state.

#### Python plots
- `plot.py [--file filename] [--elements] [--hist bins_number] [--state]` for the state vectors output file. It shows orbital elements against time, position and velocity, semi-major axis histogram at the end. 
- `display.py [--file filename] [--hist bins_width] [--div] [--r]` for the orbital elements output. It shows semi-major axis histogram, a(e) and a(w) plots. `--div` is to plot on another figure the initial distribution of the asteroids elements. `--r` shows the apoapsis(periapsis) distribution on another figure.