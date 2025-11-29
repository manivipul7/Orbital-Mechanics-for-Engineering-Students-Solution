# Space-Sciences-and-Astrodynamics

This repository contains Python and implementations of orbital mechanics and astrodynamics tools. It provides a comprehensive set of functions for simulating orbital motion, performing numerical methods, analyzing spacecraft trajectories, and modeling celestial mechanics. The repository is based on foundational principles and extends concepts found in *Orbital Mechanics for Engineering Students* by Howard D. Curtis.

## Table of Contents
- [Overview](#overview)
- [Requirements](#requirements)
- [Features](#features)
- [Course Code Modules](#course-code-modules)
- [Python Code Modules](#Python-code-modules)
- [Usage](#usage)
- [Installation](#installation)
- [Future Updates](#future-updates)

---

## Overview

The **Space Sciences and Astrodynamics** repository provides a robust toolkit for modeling, simulating, and analyzing space-based systems. It supports a wide range of functionalities including:
- Conversion between orbital elements and Cartesian state vectors.
- Simulation of spacecraft trajectories under gravitational and non-gravitational forces.
- Calculation of interplanetary trajectories and orbital perturbations.
- Numerical solutions to astrodynamics equations using advanced algorithms.
- Visualization of orbits, ground tracks, and other dynamic phenomena.

The toolkit is designed for educational and research purposes and is a resource for those interested in understanding or working within astrodynamics to explore and develop a deeper understanding of orbital mechanics.

## Requirements

### Python
- **Python**: Version 3.8 or later recommended
- Required packages:
  - `numpy`
  - `scipy`
  - `matplotlib` 


## Features
- **Coordinate Transformations**: Convert between orbital elements and state vectors.
- **Orbit Propagation**: Simulate two-body and three-body systems.
- **Trajectory Analysis**: Compute interplanetary transfers and lunar trajectories.
- **Numerical Methods**: Implement advanced solvers like Runge-Kutta, Heun’s method, and bisection.
- **Ephemeris Calculations**: Predict positions of celestial bodies for a given Julian date.
- **Perturbation Modeling**: Analyze gravitational (e.g., J2 effects) and non-gravitational forces (e.g., solar radiation pressure).

## Files and Modules

### Repository Structure
```
Space-Sciences-and-Astrodynamics/
├── .github/workflows/
├── doc/source/
├── plots/
├── curtis_scripts/
│   ├── python/
├── course_scripts/
├── tests/
│   ├── python/   
├── Algorithm_Number_Reference.md
├── CHANGELOG.md
├── CMakeLists.txt
├── LICENSE
├── pyproject.toml
├── CHANGELOG.md
└── README.md
```

### Course Code Modules

---

#### `Coordinate Transformation.py` 
- Converts between orbital elements and state vectors for orbital dynamics, calculating angular momentum, eccentricity, inclination, and other essential values.

---

#### `Family_of_orbits.py` 
- Simulates and visualizes a family of orbits with varying semi-major axis and eccentricity. Also includes an animation function for orbit visualization.

---

#### `Orbital Elements.py`
- Calculates orbital parameters from initial state vectors, such as vector magnitude, radial velocity, and angular momentum.

---

#### `Orbit_Simulator.py`
- Simulates orbital motion with transformations between orbital elements and Cartesian coordinates, providing a 3D orbit visualization.

<br>

### MATLAB Python and C++ Code Modules

---

#### `atan2d_0_360.py`

- The `atan2d_0_360` function is a specialized variant of the traditional `arctan` function, designed for scenarios where:
  - Angles need to be normalized to the range [0, 360] instead of the standard [-180, 180].
  - Consistent angle representation is required for astrodynamics, robotics, and other applications.
  - Improved readability and usability of angular calculations are beneficial in workflows.

![atan2d_0_360_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/atan2d_0_360_plot.png)

---

#### `atmosisa.py`

- The `atmosisa` function calculates atmospheric properties using an approximation of the International Standard Atmosphere (ISA) model. It provides:
  - Temperature (K) based on altitude in standard atmospheric layers.
  - Pressure (Pa) and density (kg/m³) for accurate environmental modeling.
  - Speed of sound (m/s) for applications in aerodynamics and astrodynamics.

![atmosisa_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/atmosisa_plot.png)

---

#### `atmosphere.py`

- The `atmosphere` function calculates atmospheric density for altitudes ranging from sea level to 1000 km using exponential interpolation. This provides:
  - Accurate density values based on the U.S. Standard Atmosphere 1976 model.
  - Support for high-altitude atmospheric modeling for aerospace and astrodynamics applications.
  - Smooth transitions in density across altitude intervals for improved numerical stability.

![atmosphere_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/atmosphere_plot.png)

---

#### `bisect.py`

- The `bisect` function implements the bisection method for root-finding. This numerical algorithm is designed for scenarios where:
  - The root of a continuous function needs to be located within a specified interval.
  - High precision is required, with results computed to within a user-defined tolerance.
  - The method's simplicity and robustness are ideal for solving equations in astrodynamics and related fields.

![bisect_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/bisect_plot.png)

---

#### `coe_from_sv.py`

- The `coe_from_sv` function calculates classical orbital elements (COEs) from a given state vector using Algorithm 4.1. This is essential for:
  - Translating position and velocity vectors into orbital parameters.
  - Characterizing orbits in terms of eccentricity, inclination, and other key metrics.
  - Enabling comprehensive analysis and simulation in astrodynamics.


![coe_from_sv_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/coe_from_sv_plot.png)

---

#### `cowell.py`

- The `cowell` function implements numerical integration using MATLAB's `ode45` or Python's `solve_ivp` to model satellite motion under atmospheric drag. It is designed for scenarios where:
  - Detailed trajectory predictions are needed, accounting for atmospheric drag effects.
  - Orbital evolution over extended periods must be analyzed.
  - Accurate modeling of satellite re-entry or decay due to drag is required.


![cowell_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/cowell_plot.png)

---

#### `dcm_from_q.py`

- The `dcm_from_q` function calculates the direction cosine matrix (DCM) from a quaternion. This is essential for:
  - Translating quaternion-based orientation into matrix form for applications in aerospace and robotics.
  - Enabling efficient rotation transformations in simulations and real-time systems.
  - Maintaining numerical stability and accuracy in orientation calculations.

![dcm_from_q_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/dcm_from_q_plot.png)

---

#### `dcm_to_euler.py`

- The `dcm_to_euler` function extracts the classical Euler angles from a direction cosine matrix (DCM) using the sequence $R_3(\gamma) \cdot R_1(\beta) \cdot R_3(\alpha)$. This is useful for:
  - Converting rotational transformations into intuitive angle-based representations.
  - Analyzing and visualizing orientation in astrodynamics and robotics.
  - Supporting systems requiring Euler angle decomposition from matrix representations.


![dcm_to_euler_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/dcm_to_euler_plot.png)

---

#### `dcm_to_ypr.py`

- The `dcm_to_ypr` function extracts the yaw, pitch, and roll angles from a direction cosine matrix (DCM) using the sequence $$R_1(\text{roll}) \cdot R_2(\text{pitch}) \cdot R_3(\text{yaw}) $$. This is useful for:
  - Converting matrix-based rotations into intuitive yaw, pitch, and roll angles.
  - Applications in aerospace, robotics, and navigation systems that require angle-based orientation.
  - Analyzing orientation data in a more comprehensible format.

![dcm_to_ypr_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/dcm_to_ypr_plot.png)

---

#### `encke.py`

- The `encke` function uses Encke's method to numerically integrate the equations of motion for a satellite under the influence of $J_2$ gravitational perturbations. This code is designed for scenarios where:
  - Precise orbit propagation is required over long periods while accounting for perturbative forces.
  - Modeling the effects of Earth's oblateness $$( J_2 )$$ on orbital elements is essential.
  - Computational efficiency is important for handling small perturbations in satellite motion.


![encke_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/encke_plot1.png)

![encke_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/encke_plot2.png)

![encke_plot3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/encke_plot3.png)

---

#### `f_and_g.py`

- The `f_and_g` function calculates the Lagrange coefficients $f$ and $g$, which are essential for solving orbital propagation problems. These coefficients:
  - Aid in determining the position and velocity of a celestial body after a time interval.
  - Incorporate the universal anomaly $x$ and other orbital parameters for precise calculations.
  - Support astrodynamics workflows requiring accurate orbital state predictions.

![f_and_g_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/f_and_g_plot.png)

---

#### `f_and_g_ta.py`

- The `f_and_g_ta` function calculates the Lagrange $f$ and $g$ coefficients based on the change in true anomaly ($\Delta \theta$). This function is particularly useful for:
  - Propagating orbits using true anomaly changes without requiring a full numerical integration.
  - Efficiently determining position and velocity at a new orbital location.
  - Supporting analytical orbital mechanics workflows.


![f_and_g_ta_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/f_and_g_ta_plot.png)

---

#### `fDot_and_gDot.py`

- The `fDot_and_gDot` function calculates the time derivatives of the Lagrange $f$ and $g$ coefficients, which are critical for:
  - Determining changes in position and velocity over time in orbital mechanics.
  - Supporting precise orbit propagation and trajectory predictions.
  - Analyzing dynamic changes in orbital elements using universal anomaly and radial positions.

![fDot_and_gDot_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/fDot_and_gDot_plot.png)

---

#### `fDot_and_gDot_ta.py`

- The `fDot_and_gDot_ta` function calculates the time derivatives of the Lagrange $f$ and $g$ coefficients based on a change in the true anomaly ($\Delta \theta$). This is particularly useful for:
  - Analyzing time-dependent changes in orbital position and velocity.
  - Supporting orbital mechanics calculations involving angular momentum and radial velocity.
  - Providing efficient computations for trajectory propagation in astrodynamics.

![fDot_and_gDot_ta_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/fDot_and_gDot_ta_plot.png)

---

#### `gauss.py`

- The `gauss` function implements Gauss's method for preliminary orbit determination, including iterative improvement. This algorithm calculates the state vector (position and velocity) of an orbiting body from angles-only observations at three closely spaced times. Key features include:
  - Precise orbit determination using minimal observation data.
  - Iterative refinement to improve the accuracy of the calculated state vector.
  - Support for real-world applications in satellite tracking and orbital mechanics.

![gauss_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/gauss_plot.png)

---

#### `gibbs.py`

- The `gibbs` function implements Gibbs's method for preliminary orbit determination using three coplanar geocentric position vectors. It is designed for:
  - Calculating the velocity corresponding to the second position vector.
  - Determining orbital elements such as angular momentum, eccentricity, and inclination.
  - Supporting applications in satellite tracking and celestial mechanics.

![gibbs_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/gibbs_plot.png)

---

#### `gravity_turn.py`

- The `gravity_turn` function numerically integrates the equations of motion for a gravity turn trajectory. This is particularly useful for:
  - Simulating launch vehicle trajectories under the influence of thrust, drag, and gravity.
  - Calculating dynamic parameters such as altitude, velocity, and flight path angle.
  - Analyzing velocity losses due to drag and gravity.

![gravity_turn_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/gravity_turn_plot1.png)

![gravity_turn_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/gravity_turn_plot2.png)

---

#### `ground_track.py`

- The `ground_track` function calculates and plots the ground track of a satellite using its orbital elements. This is useful for:
  - Visualizing the satellite's trajectory relative to Earth's surface over multiple orbits.
  - Simulating orbital evolution, including nodal regression and perigee precession due to $J_2$ perturbations.
  - Analyzing the satellite's coverage and ground visibility.

![ground_track_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/ground_track_plot.png)

---

#### `heun.py`

- The `heun` function implements Heun's method, a predictor-corrector numerical integration technique, for solving systems of first-order ordinary differential equations (ODEs). This function is designed for:
  - Accurate numerical integration of time-dependent ODEs in engineering and physics applications.
  - Handling systems with stringent error tolerance and iterative correction.
  - Supporting user-defined differential equations through callable functions.

![heun_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/heun_plot.png)

---

#### `integrate_thrust.py`

- The `integrate_thrust` function uses numerical integration to model the dynamics of a spacecraft during a delta-v burn. It calculates the resulting trajectory and orbital parameters after the burn. This function is ideal for:
  - Simulating orbital maneuvers, including changes in velocity and trajectory.
  - Computing the apogee and other orbital elements of the post-burn orbit.
  - Evaluating thrust, mass loss, and other dynamic factors affecting spacecraft motion.

![integrate_thrust_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/integrate_thrust_plot.png)

---

#### `interplanetary.py`

- The `interplanetary` function calculates the trajectory of a spacecraft traveling from one planet to another using patched conic approximations. This function is ideal for:
  - Planning interplanetary missions by determining heliocentric state vectors at departure and arrival.
  - Computing the spacecraft's hyperbolic excess velocities at both planets.
  - Supporting mission design and analysis with key orbital elements and time-of-flight calculations.

![interplanetary_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/interplanetary_plot.png)

---

#### `J0.py`

- The `J0` function calculates the Julian day number at 0 UT (Universal Time) for a given date. This is essential for:
  - Converting calendar dates into a numerical format suitable for astronomical and orbital calculations.
  - Supporting time-dependent computations in astrodynamics and celestial mechanics.
  - Ensuring compatibility with standard equations for Julian date arithmetic.

![J0_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/J0_plot.png)

---

#### `J2_perturbation.py`

- The `J2_perturbation` function numerically integrates Gauss's planetary equations to analyze the effects of Earth's oblateness ($J_2 $) on the orbital elements of a satellite. This function is designed for:
  - Investigating long-term perturbations in orbital parameters caused by $J_2$.
  - Simulating changes in right ascension, inclination, eccentricity, and argument of perigee.
  - Supporting mission planning and satellite lifetime analysis with precise perturbation modeling.

![J2_perturbation_plot_1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/J2_perturbation_plot_1.png)

![J2_perturbation_plot_2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/J2_perturbation_plot_2.png)

![J2_perturbation_plot_3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/J2_perturbation_plot_3.png)

---

#### `kepler_E.py`

- The `kepler_E` function solves Kepler's equation $E - e \sin(E) = M$ for the eccentric anomaly ($E$) using Newton's method. This function is essential for:
  - Determining orbital positions in elliptical trajectories.
  - Supporting astrodynamics calculations where precise orbital element computation is required.
  - Converting mean anomaly ($M$) to eccentric anomaly ($E$) for orbital mechanics analyses.

![kepler_E_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/kepler_E_plot.png)

---

#### `kepler_H.py`

- The `kepler_H` function solves Kepler's equation for hyperbolic orbits, $e \sinh(F) - F = M$, using Newton's method to compute the hyperbolic eccentric anomaly ($F$). This function is essential for:
  - Determining orbital positions in hyperbolic trajectories.
  - Supporting interplanetary mission design involving escape or flyby trajectories.
  - Converting hyperbolic mean anomaly ($M$) to hyperbolic eccentric anomaly ($F$).

![kepler_H_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/kepler_H_plot.png)

---

#### `kepler_U.py`

- The `kepler_U` function solves the universal Kepler's equation to compute the universal anomaly ($x$) using Newton's method. This function is crucial for:
  - Propagating orbits in any conic section (elliptical, parabolic, or hyperbolic).
  - Calculating position and velocity vectors over time for celestial bodies or spacecraft.
  - Enabling generalized solutions for orbital mechanics problems independent of orbit type.

![kepler_U_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/kepler_U_plot.png)

---

#### `lambert.py`

- The `lambert` function solves Lambert's problem, which calculates the velocity vectors required to transfer a spacecraft between two position vectors in a specified time. This is essential for:
  - Planning interplanetary and orbital transfer trajectories.
  - Computing initial and final velocities for given orbital configurations.
  - Supporting mission design and analysis in astrodynamics.

![lambert_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/lambert_plot.png)

---

#### `los.py`

- The `los` function determines whether a satellite is in Earth's shadow based on the satellite's and Sun's Earth-Centered Inertial (ECI) position vectors. This function is critical for:
  - Assessing solar exposure for satellites in orbit.
  - Determining periods of Earth shadow for satellite operations.
  - Supporting power and thermal analysis for spacecraft.

![los_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/los_plot.png)

---

#### `LST.py`

- The `LST` function calculates the Local Sidereal Time (LST) at a specified location on Earth. This is essential for:
  - Determining the orientation of the Earth relative to celestial objects.
  - Supporting astronomical observations and satellite tracking.
  - Converting time and geographic location into a sidereal time format for astrodynamics computations.

![LST_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/LST_plot.png)

---

#### `lunar_perturbation_test.py`

- The `lunar_perturbation_test` function models the effects of the Moon's gravitational perturbations on satellite orbits using the Gauss variational equations. This function is ideal for:
  - Analyzing long-term orbital changes due to lunar gravity.
  - Simulating orbital element variations such as inclination, right ascension, and argument of perigee.
  - Supporting mission planning and stability assessments for Earth-orbiting satellites.

![lunar_perturbation_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/lunar_perturbation_plot1.png)

![lunar_perturbation_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/lunar_perturbation_plot2.png)

---

#### `lunar_position.py`

- The `lunar_position` function calculates the geocentric equatorial position vector of the Moon for a given Julian date. This is essential for:
  - Determining the Moon's position relative to Earth in celestial mechanics.
  - Supporting lunar mission planning and astrodynamics calculations.
  - Providing precise Moon position data for time-sensitive orbital computations.

![lunar_position_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/lunar_position_plot.png)

---

#### `lunar_trajectory.py`

- The `lunar_trajectory` function models the trajectory of a spacecraft under the gravitational influence of Earth and the Moon using a restricted three-body problem framework. Key features include:
  - Accurate modeling of spacecraft motion with initial conditions and parameters.
  - Visualization of the trajectory in both the Earth-Centered Inertial (ECI) frame and a Moon-fixed frame.
  - Calculation of key points, such as perilune and lunar arrival conditions, for trajectory analysis.
  - Support for mission planning, trajectory optimization, and astrodynamics research.

![lunar_trajectory_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/lunar_trajectory_plot1.png)

![lunar_trajectory_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/lunar_trajectory_plot2.png)

![lunar_trajectory_plot3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/lunar_trajectory_plot3.png)

---

#### `month_planet_names.py`

- The `month_planet_names` function returns the name of a month and a planet based on their respective numeric IDs. This utility is designed for:
  - Mapping numeric IDs to human-readable month and planet names.
  - Supporting educational tools or software requiring quick access to such mappings.
  - Simplifying conversions between numerical representations and their corresponding names.

![month_planet_names_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/month_planet_names_plot.png)

---

#### `orbit.py`

- The `orbit` function computes the orbit of a spacecraft by numerically solving the two-body problem using the Runge-Kutta-Fehlberg (RKF45) method. This function is designed for:
  - Simulating orbital trajectories and calculating spacecraft positions and velocities over time.
  - Determining key orbital parameters, including maximum and minimum radii and the corresponding velocities.

![orbit_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/orbit_plot.png)

---

#### `planet_elements_and_sv.py`

- The `planet_elements_and_sv` function calculates the heliocentric orbital elements and state vectors (position and velocity) of a planet for a given date and time. This function is essential for:
  - Determining the current position and motion of planets in the solar system.
  - Supporting astrodynamics and celestial mechanics computations.
  - Providing data for mission planning and planetary observation.

![planet_elements_and_sv_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/planet_elements_and_sv_plot.png)

---

#### `q_from_dcm.py`

- The `q_from_dcm` function calculates the quaternion representation of a given direction cosine matrix (DCM). This is essential for:
  - Converting rotation matrices into quaternions for use in spacecraft attitude dynamics and robotics.
  - Enabling efficient computation of orientation changes.
  - Supporting applications in control systems, simulation, and 3D transformations.

![q_from_dcm_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/q_from_dcm_plot.png)

---

#### `quat_rotate.py`

- The `quat_rotate` function rotates a vector by a unit quaternion, providing a method for:
  - Transforming vectors in 3D space using quaternion-based rotation.
  - Supporting spacecraft attitude simulations and 3D graphics applications.
  - Simplifying rotational transformations with efficient quaternion mathematics.

#### `quatinv.py` and `quatmultiply.py`

- These Python scripts implement the functionality of MATLAB's built-in `quatinv` and `quatmultiply` functions. They are essential for:
  - Calculating the inverse of a quaternion with `quatinv`.
  - Performing quaternion multiplication with `quatmultiply`.

- These scripts replicate MATLAB's behavior, enabling seamless transition for MATLAB users moving to Python environments. They are utilized in functions such as `quat_rotate` for quaternion-based calculations.

![quat_rotate_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/quat_rotate_plot.png)

---

#### `ra_and_dec_from_r.py`

- The `ra_and_dec_from_r` function calculates the right ascension and declination from a geocentric equatorial position vector. This function is essential for:
  - Converting position vectors into angular coordinates for use in celestial mechanics and navigation.
  - Determining the orientation of an object in Earth's sky from its position in space.
  - Supporting astrodynamics applications, including satellite tracking and orbital visualization.

![ra_and_dec_from_r_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/ra_and_dec_from_r_plot.png)

---

#### `relative_motion.py`

- The `relative_motion` function calculates and plots the motion of a chaser satellite (B) relative to a target satellite (A) in orbit. This function is useful for:
  - Simulating relative motion dynamics between two satellites for formation flying or rendezvous missions.
  - Visualizing relative trajectories in a co-moving reference frame.
  - Supporting analysis of proximity operations in astrodynamics.

![relative_motion_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/relative_motion_plot.png)

---

#### `rk1_4.py`

- The `rk1_4` function implements the Runge-Kutta numerical integration methods (RK1, RK2, RK3, and RK4) to solve systems of first-order ordinary differential equations $\frac{dy}{dt} = f(t, y)$. This function is essential for:
  - Performing numerical simulations of dynamic systems.
  - Providing flexibility to choose between integration methods depending on accuracy and computational efficiency.
  - Supporting the analysis of time-dependent problems in astrodynamics and engineering.

![rk1_4_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/rk1_4_plot.png)

---

#### `rkf45.py`

- The `rkf45` function implements the Runge-Kutta-Fehlberg 4(5) method, a numerical integration algorithm with adaptive step size control. It is designed for:
  - Solving systems of first-order ordinary differential equations $\frac{dy}{dt} = f(t, y)$ with high accuracy and efficiency.
  - Applications in astrodynamics, engineering, and physics requiring precise integration over varying time steps.
  - Supporting problems with stringent error tolerance and dynamic step size adjustments.

![rkf45_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/rkf45_plot.png)

![rkf45_lunar_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/rkf45_lunar_plot.png)

---

#### `rv_from_observe.py`

- The `rv_from_observe` function calculates the geocentric equatorial position and velocity vectors of an object from radar observations, including range, azimuth, elevation angle, and their rates. This function is essential for:
  - Translating radar observation data into state vectors for orbit determination.
  - Supporting satellite tracking and orbit prediction applications.
  - Providing a foundation for converting observational data into orbital elements.

![rv_from_observe_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/rv_from_observe_plot.png)

---

#### `rv_from_r0v0.py`

- The `rv_from_r0v0` function computes the state vector (position and velocity) of a spacecraft at a future time, given its initial state vector and the elapsed time. This is fundamental for:
  - Propagating orbital trajectories using Kepler's equations and universal anomaly methods.
  - Simulating spacecraft motion in the two-body problem framework.
  - Supporting orbital dynamics analysis and mission planning tasks.


![rv_from_r0v0_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/rv_from_r0v0_plot.png)

---

#### `rv_from_r0v0_ta.py`

- The `rv_from_r0v0_ta` function calculates the state vector (position and velocity) of a spacecraft after a specified change in true anomaly, given the initial state vector. This function is essential for:
  - Propagating orbital trajectories based on changes in true anomaly.
  - Determining spacecraft position and velocity at specific points in the orbit.
  - Supporting astrodynamics analysis and mission design tasks.

![rv_from_r0v0_ta_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/rv_from_r0v0_ta_plot.png)

---

#### `rva_relative.py`

- The `rva_relative` function calculates the position, velocity, and acceleration of a spacecraft B relative to spacecraft A in A's Local Vertical Local Horizontal (LVLH) frame. This function is vital for:
  - Analyzing relative motion in rendezvous and proximity operations.
  - Simulating formation flying and cooperative satellite missions.
  - Supporting mission planning for docking and station-keeping tasks.

![rva_relative_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/rva_relative_plot.png)

---

#### `simpsons_lunar_ephemeris.py`

- The `simpsons_lunar_ephemeris` function computes the state vector (position and velocity) of the Moon relative to Earth's geocentric equatorial frame for a given Julian date. This function is based on a curve fit to JPL's DE200 ephemeris model, providing:
  - Accurate lunar position and velocity data for astrodynamics applications.
  - Essential inputs for spacecraft mission planning and lunar navigation.
  - A simplified, efficient approach suitable for onboard flight software.

![simpsons_lunar_ephemeris_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/simpsons_lunar_ephemeris_plot1.png)

---

#### `solar_perturbation.py`

- The `solar_perturbation` function simulates the effects of solar gravitational perturbation on a spacecraft's orbit by integrating the Gauss variational equations. This function is designed for:
  - Analyzing the long-term influence of solar gravity on satellite orbits.
  - Supporting mission planning by quantifying perturbations on different orbital regimes (e.g., LEO, HEO, GEO).
  - Providing insights into orbital element variations due to external gravitational effects.

![solar_perturbation_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/solar_perturbation_plot1.png)

![solar_perturbation_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/solar_perturbation_plot2.png)

---

#### `solar_position.py`

- The `solar_position` function calculates the geocentric equatorial position vector of the Sun for a given Julian date. This is crucial for:
  - Determining the Sun's position relative to Earth in celestial mechanics and satellite mission planning.
  - Supporting solar radiation pressure calculations and perturbation modeling.
  - Providing accurate inputs for astrodynamics and space mission design.

![solar_position_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/solar_position_plot.png)

---

#### `solar_radiation_pressure.py`

- The `solar_radiation_pressure` function models the effects of solar radiation pressure on a satellite's orbital elements by solving the Gauss planetary equations. This is particularly useful for:
  - Evaluating the long-term influence of solar radiation on satellite orbits.
  - Supporting mission planning by quantifying non-gravitational perturbations.
  - Analyzing changes in angular momentum, eccentricity, and other orbital parameters over time.

![solar_radiation_pressure_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/solar_radiation_pressure_plot.png)

![solar_radiation_pressure_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/solar_radiation_pressure_plot2.png)

---

#### `spinning_top.py`

- The `spinning_top` function numerically integrates Euler's equations of motion for a spinning top. It uses quaternions to track the time evolution of the top's orientation. This implementation is designed for:
  - Understanding the rotational dynamics of a rigid body under gravity.
  - Analyzing Euler angles, angular velocities, and the effects of moments of inertia.
  - Providing insights into the nutation, precession, and spin behavior of spinning tops.

![spinning_top_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/spinning_top_plot1.png)

![spinning_top_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/spinning_top_plot2.png)

![spinning_top_plot3.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/spinning_top_plot3.png)

---

#### `stumpC.py`

- The `stumpC` function evaluates the Stumpff function $C(z)$, a key mathematical function in orbital mechanics, particularly for solving Kepler's equation using the universal variable formulation. This implementation is designed for:
  - Efficiently handling cases for positive, negative, and zero values of $z$.
  - Supporting numerical algorithms in astrodynamics that require Stumpff function computations.
  - Providing robust and reliable results for diverse orbital scenarios.

![stumpC_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/stumpC_plot.png)

---

#### `stumpS.py`

- The `stumpS` function computes the Stumpff function $S(z)$, essential in solving Kepler's equation using the universal variable formulation. It handles different cases of $z$ (positive, negative, and zero) to provide:
  - Accurate evaluations for orbital mechanics problems.
  - Support for numerical algorithms in astrodynamics requiring Stumpff function computations.
  - Reliable results across various orbital scenarios.

![stumpS_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/stumpS_plot.png)

---

#### `sv_from_coe.py`

- The `sv_from_coe` function computes the state vector (position and velocity) from classical orbital elements. This is particularly useful in astrodynamics for:
  - Translating orbital elements like eccentricity, inclination, and true anomaly into precise position and velocity vectors in a geocentric equatorial frame.
  - Supporting simulations and analyses of orbital mechanics, mission planning, and trajectory design.
  - Ensuring consistent and accurate transformations between orbital and Cartesian coordinates.

![sv_from_coe_plot.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/sv_from_coe_plot.png)

---

#### `threebody.py`

- The `threebody` function provides a graphical simulation of the motion of three celestial bodies interacting under gravitational forces in a two-dimensional plane. It uses:
  - Newton's laws of motion and the gravitational constant to compute positions and velocities over time.
  - Numerical integration (Runge-Kutta in MATLAB and Euler's method in Python) to solve the system of differential equations governing the bodies' motion.

![threebody_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/threebody_plot1.png)

![threebody_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/threebody_plot2.png)

---

#### `twobody3d.py`

- The `twobody3d` function numerically solves the two-body problem in three dimensions relative to an inertial frame. The simulation uses the RKF4(5) numerical integration method to model the gravitational interactions between two massive bodies.
  - Inputs include initial positions and velocities of the bodies, their masses, and the time interval of the simulation.
  - Outputs provide trajectories, velocities, and the center of mass for the system at discrete time intervals.

![twobody3d_plot1.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/twobody3d_plot1.png)

![twobody3d_plot2.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/twobody3d_plot2.png)

![twobody3d_plo3t.png](https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics/blob/main/plots/twobody3d_plot3.png)

---

## Usage
- Build Python: `pip install .`

### Python
Clone the repository and navigate to the python files directory. Install the required .py packages and run scripts directly.

## Installation

Clone the repository:
```bash
git clone https://github.com/Keeby-Astro/Space-Sciences-and-Astrodynamics.git
```

---