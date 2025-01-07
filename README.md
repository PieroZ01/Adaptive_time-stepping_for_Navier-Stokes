<a name="readme-top"></a>

[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
[![Gmail][gmail-shield]][gmail-url]

<br />
<div align="center">

<h1 align="center">Adaptive time-stepping for Navier-Stokes</h1>
<h3 align="center">Advanced Topics in Scientific Computing Course Final Project</h3>
<h3 align="center">SDIC Master Degree, University of Trieste (UniTS) & SISSA</h3>
<h3 align="center">2024-2025</h3>

<p align="center">
    Adaptive time-stepping for solving the incompressible Navier-Stokes equations with a finite element method using the <a href="https://www.dealii.org">deal.II</a> library.
    <br />
    <br />
</div>

<!-- TABLE OF CONTENTS -->
<div style="width: 360px; text-align: center; border: 2px solid currentColor; padding: 10px 10px 10px 10px; border-radius: 10px; margin: auto;">
  <h4>📑 Table of Contents</h4>
  <ul style="list-style-type: none; padding: 0;">
    <li><a href="#author-info">Author Info</a></li>
    <li><a href="#project-overview">Project Overview</a></li>
    <li><a href="#adaptive-time-stepping-algorithm">Adaptive Time-Stepping Algorithm</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#results">Results</a></li>
    <li><a href="#conclusions">Conclusions</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#references">References</a></li>
  </ul>
</div>
</br>

<!-- AUTHORS INFO-->
## Author Info

| Name | Surname | Student ID | UniTs email | Personal email | Master course |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Piero | Zappi | SM3600004 | `piero.zappi@studenti.units.it` | `piero.z.2001@gmail.com` | SDIC |

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Project Overview

### Project Goal

The aim of this project is to implement an **adaptive time-stepping algorithm** for solving the incompressible Navier-Stokes equations with a finite element method using the [deal.II](https://www.dealii.org) library. The project is developed as the final project for the *Advanced Topics in Scientific Computing* course of the *<a href="https://sdic.units.it">SDIC</a>* Master Degree at the *University of Trieste (UniTS)* and *SISSA* for the academic year 2024-2025.

### Introduction

**Adaptive time-stepping** is an important technique in the numerical solution of time-dependent partial differential equations (PDEs) for controlling the **accuracy** of the simulations and reducing the **computational costs**. It is particularly useful for problems with a wide range of time scales, where a fixed time step size would be inefficient or inaccurate.
An adaptive time-stepping algorithm adjusts the time step size at each time step of the simulation based on some error estimator, which suggests how much the time step should be increased or decreased to achieve a desired level of accuracy: basically, the time step is increased when the error is small and decreased when the error is large. The dynamic adjustment of the time step size allows the simulation to be more efficient and accurate; in fact, if the time step size is too small and the solution is already accurate, a lot of unnecessary computations are performed, resulting in a waste of computational resources, making the simulation inefficient. On the other hand, if the time step size is too large, the solution may become inaccurate. Therefore, an adaptive time-stepping algorithm aims to find the optimal time step size that minimizes the computational cost of a simulation, while maintaining the desired level of accuracy: this is achieved by dynamically adjusting the time step size based on the error estimator during the simulation, following the time evolution of the physics of the given problem.

### Project Foundation

The development of the project is based on the program of the **Step-35 tutorial** of the deal.II library [<a href="#ref1">1</a>], which solves the incompressible time-dependent Navier-Stokes equations using a projection method, which decouples the velocity and pressure fields, introducing auxiliary variables. For details on the implementation of the projection method, the fully discrete setting and the code structure, please refer to the [Step-35 tutorial](https://dealii.org/current/doxygen/deal.II/step_35.html).
This program provides a solid foundation for solving the Navier-Stokes equations and can be extended to include adaptive time-stepping capabilities. In fact, the original program uses a fixed time step size, which can be replaced by an adaptive time-stepping algorithm to improve the efficiency and the accuracy of the simulations. The code has then been modified and extended to work with a variable time step size and to include an adaptive time-stepping technique, which dinamically adjusts the time step size based on an error estimator, following the time evolution of the solution of the Navier-Stokes equations during the simulation.


## Adaptive Time-Stepping Algorithm

## Usage

## Results

## Conclusions


<!-- LICENSE -->
## License

Distributed under the MIT License. See [`LICENSE`](./LICENSE) for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- REFERENCES -->
## References

<a id="ref1"></a>
[1] <a href="https://dealii.org/current/doxygen/deal.II/step_35.html">Step-35 tutorial program</a> of the deal.II library: a projection solver for the Navier–Stokes equations.

<a id="ref2"></a>
[2] Berrone, S. and Marro, M., 2009. Space–time adaptive simulations for unsteady Navier–Stokes problems. Computers & Fluids, 38(6), pp.1132-1144.

<a id="ref3"></a>
[3] Kay, D.A., Gresho, P.M., Griffiths, D.F. and Silvester, D.J., 2010. Adaptive time-stepping for incompressible flow part II: Navier–Stokes equations. SIAM Journal on Scientific Computing, 32(1), pp.111-128.

<a id="ref4"></a>
[4] Boisneault, A., Dubuis, S. and Picasso, M., 2023. An adaptive space-time algorithm for the incompressible Navier-Stokes equations. Journal of Computational Physics, 493, p.112457.

<a id="ref5"></a>
[5] John, V. and Rang, J., 2010. Adaptive time step control for the incompressible Navier–Stokes equations. Computer Methods in Applied Mechanics and Engineering, 199(9-12), pp.514-524.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- Contacts -->
[license-shield]: https://img.shields.io/badge/License-MIT-blue.svg?style=for-the-badge
[license-url]: https://github.com/PieroZ01/Adaptive_time-stepping_for_Navier-Stokes/blob/master/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-blue?style=for-the-badge&logo=linkedin&logoColor=white&colorB=0077B5
[linkedin-url]: https://www.linkedin.com/in/pierozappi/
[gmail-shield]: https://img.shields.io/badge/-Gmail-red?style=for-the-badge&logo=gmail&logoColor=white&colorB=red
[gmail-url]: mailto:piero.z.2001@gmail.com