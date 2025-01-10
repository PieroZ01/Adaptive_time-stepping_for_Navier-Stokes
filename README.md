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
    Adaptive time-stepping for solving the incompressible Navier-Stokes equations with a finite element method using the <a href="https://www.dealii.org">deal.II</a> library
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

| Name | Surname | Student ID | UniTS email | Personal email | Master course |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Piero | Zappi | SM3600004 | `piero.zappi@studenti.units.it` | `piero.z.2001@gmail.com` | SDIC |

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Project Overview

### Project Objective

The aim of this project is to implement an **adaptive time-stepping algorithm** for solving the incompressible Navier-Stokes equations with a finite element method using the [deal.II](https://www.dealii.org) library. The project is developed as the final project for the *Advanced Topics in Scientific Computing* course of the *Scientific and Data-Intensive Computing* Master Degree at the *University of Trieste (UniTS)* and *SISSA* for the academic year 2024-2025.

### Introduction

**Adaptive time-stepping** is an important technique in the numerical solution of time-dependent partial differential equations (PDEs) for controlling the **accuracy** of the simulations and reducing the **computational costs**. It is particularly useful for problems with a wide range of time scales, where a fixed time step size would be inefficient or inaccurate.
An adaptive time-stepping algorithm adjusts the time step size at each time step of the simulation based on some error estimator, which suggests how much the time step should be increased or decreased to achieve a desired level of accuracy: basically, the time step is increased when the error is small and decreased when the error is large. The dynamic adjustment of the time step size allows the simulation to be more efficient and accurate; in fact, if the time step size is too small and the solution is already accurate, a lot of unnecessary computations are performed, resulting in a waste of computational resources, making the simulation inefficient. On the other hand, if the time step size is too large, the solution may become inaccurate. Therefore, an adaptive time-stepping algorithm aims to find the optimal time step size that minimizes the computational cost of a simulation, while maintaining the desired level of accuracy: this is achieved by dynamically adjusting the time step size based on the error estimator during the simulation, following the time evolution of the physics of the given problem.

### Project Foundation

The development of the project is based on the program of the **Step-35 tutorial** of the deal.II library [<a href="#ref1">1</a>], which solves the incompressible time-dependent Navier-Stokes equations using a **projection method**, which decouples the velocity and pressure fields, introducing auxiliary variables. For details on the implementation of the projection method, the fully discrete setting and the code structure, please refer to the [Step-35 tutorial](https://dealii.org/current/doxygen/deal.II/step_35.html).
This program provides a solid foundation for solving the Navier-Stokes equations and can be extended to include adaptive time-stepping capabilities. In fact, the original program uses a fixed time step size, which can be replaced by an adaptive time-stepping algorithm to improve the efficiency and the accuracy of the simulations. The code has then been modified and extended to work with a variable time step size and to include an adaptive time-stepping technique, which dinamically adjusts the time step size based on an error estimator, following the time evolution of the solution of the Navier-Stokes equations during the simulation.

### Project Structure

The project is structured as follows:

```bash
.
├── CMakeLists.txt # CMake configuration file for compiling the project
├── LICENSE # MIT License file
├── README.md # This README file
├── adaptive_NS.cc # Main program file for solving the Navier-Stokes equations with adaptive time-stepping
├── data # Folder containing the results of some simulations
│   ├── Re100_norm_level1
│   ├── Re100_norm_level2
│   ├── Re200_norm_level2
│   └── Re50_norm_level2
├── images # Folder containing some images
│   ├── sol.gif
│   └── sol_LIC.gif
├── nsbench.inp # Input file containing the geometry of the test case
├── parameter-file.prm # Input file containing the parameters for the simulations
├── time-step_analysis.ipynb # Jupyter notebook for analyzing the time step size evolution during the simulations
└── time_steps.csv # CSV file for saving the time step size evolution during a simulation
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Adaptive Time-Stepping Algorithm

As stated in the Step-35 tutorial of the deal.II library, currently there is not much theory about projection methods with variable time steps. However, in literature there are some works that have implemented adaptive time-stepping algorithms for the Navier-Stokes equations, such as the works of Berrone and Marro [<a href="#ref2">2</a>], Kay et al. [<a href="#ref3">3</a>], Boisneault et al. [<a href="#ref4">4</a>] and John and Rang [<a href="#ref5">5</a>]. These works provide a solid foundation for the implementation of adaptive time-stepping algorithms for the Navier-Stokes equations, which can be used as a reference for the development of an adaptive time-stepping algorithm in this project.

Mainly, there are **two different approaches** for implementing adaptive time-stepping algorithms for the Navier-Stokes equations. The first one develops an adaptive time step control based on **two solutions** which are computed with schemes of **different orders** of accuracy. The **difference** between these two solutions can be used to estimate an appropriate time step size for the next time step.

### First Approach

We now try to summarize and schematize the way in which this first type of adaptive time-stepping algorithm works: the difference $e$ between the two solutions is measured with some functional, such as the $L^2$ norm. Assuming that the less accurate solution $u$ has an error of order $\mathcal{O}(\Delta t^p)$ and the more accurate solution $\bar{u}$ has an error of order $\mathcal{O}(\Delta t^{p+1})$, then the error $e$ is of order $\mathcal{O}(\Delta t^p)$:

$$
e = ||\bar{u}-u||_{L^2} \approx C \Delta t^p
$$

Given a desired tolerance $\epsilon>0$, the following relation can be used to estimate the optimal time step size $\Delta t_{\text{opt}}$ for the next time step:

$$
\frac{e}{\epsilon} = \frac{\Delta t^p}{\Delta t_{\text{opt}}^p}
$$

The optimal time step size $\Delta t_{\text{opt}}$ for the next time step is then computed as:

$$
\Delta t_{\text{opt}} = \rho \Delta t \left( \frac{\epsilon}{e} \right)^{1/p}
$$

where $\rho\in(0,1]$ is a safety factor. The optimal time step size $\Delta t_{\text{opt}}$ is then used for the next time step of the simulation.

The main drawback of this approach is its high **computational cost**, as it requires the computation of two solutions with different schemes, which can be computationally expensive. In fact, the second scheme is used only to determine the size of the next time step, which can be considered a waste of computational resources, as the cost of the simulation is practically doubled. Therefore, this approach is not very efficient and risks to nullify the benefits in terms of computational efficiency deriving from the use of an adaptive time-stepping algorithm. However, there exist some methods that can be used to compute the two solutions with a reduced computational cost solving this issue (see the works of Kay et al. [<a href="#ref3">3</a>] and John and Rang [<a href="#ref5">5</a>] for more details).

### Second Approach (Implemented in this Project)

In theory the second approach, which is actually the one *implemented in this project*, is more efficient than the first one, as it uses only **one solution** (in this specific case, the one provided by the projection method) to estimate the optimal time step size for the next time step. Basically, the time step is chosen on the basis of comparing the **change in the solution of two consecutive time steps** in the $L^2$ norm of the space-time interval. The time step size is adjusted based on the relative change in the solution, measured with a *time error indicator*, which is used to determine the optimal time step size for the next time step.

This approach, discussed in the works of Berrone and Marro [<a href="#ref2">2</a>] and Boisneault et al. [<a href="#ref4">4</a>], enables the time step size to be adjusted dynamically during the simulation based on the time evolution of the system, *"following" the physics of the problem*. The time error indicator is used to modify the time step size accordingly, ensuring that the simulation is accurate and efficient. The time step size is reduced when the system is evolving rapidly and the solution is changing significantly (i.e., the error is large), while it is increased when the system is evolving slowly and the solution is changing less (i.e., the error is small).

The following quantity has been used as a time error indicator in this project:

$$
\begin{equation*}
    (\eta^n)^2 = \sum_{K\in T_h} \left( \int_{t^{n-1}}^{t^n} \left(||\nabla\left(u_h^{n}-u_h^{n-1}\right)||^2_{L^2(K)}\right)dt \right)
\end{equation*}
$$

where $T_h$ is the mesh, $u_h^n$ is the solution at time $t^n$ and $u_h^{n-1}$ is the solution at time $t^{n-1}$. The time error indicator $\eta^n$ allows to evaluate the quality of the time step size $\Delta t^n$ according to the following principles: if $\eta^n$ is large, the time step size $\Delta t^n$ is too large and should be reduced; otherwise if $\eta^n$ is small, the time step size $\Delta t^n$ is too small and should be increased.

The following term has been introduced to normalize the time error indicator:

$$
(\sigma^n)^2 = \int_{t^{n-1}}^{t^n} |u_h^{n}|^2_{1}dt
$$

where $|u_h^{n}|^2_{1}$ is the $H^1$ semi-norm of the solution $u_h^n$.

In order to obtain at each time step the optimal time step size $\Delta t_{\text{opt}}$ for the next time step, the adaptive time-stepping algorithm employs the following bound:

$$
\begin{equation}
    (1-\alpha)\text{TOL}\leq\frac{\eta^n}{\sigma^n}\leq(1+\alpha)\text{TOL}
\end{equation}
$$

where $\text{TOL}>0$ is a given tolerance and $\alpha$ is a parameter in the range $(0,1]$. Specifically, the tolerance was set to $\text{TOL}=0.1$ and the parameter $\alpha$ was set to $\alpha=0.5$ in this project.

If the double inequality (1) is satisfied, the current time step size $\Delta t^n$ is accepted and is kept for the next time step. Otherwise, the time step size is modified according to the following rule: it is reduced if the time error indicator exceeds the upper bound and it is increased if it is below the lower bound. In particular, the optimal time step size $\Delta t_{\text{opt}}$ for the next time step is computed as:

$$
\begin{equation*}
\begin{aligned}
\text{if } (\eta^n)^2 &> (\sigma^n)^2(1+\alpha)^2\text{TOL}^2 \\
\quad \Delta t_{\text{opt}} &= \frac{\Delta t^n}{\rho_1} \\
\text{else if } (\eta^n)^2 &< (\sigma^n)^2(1-\alpha)^2\text{TOL}^2 \\
\quad \Delta t_{\text{opt}} &= \frac{\Delta t^n}{\rho_2} \\
\text{else} \\
\quad \Delta t_{\text{opt}} &= \Delta t^n
\end{aligned}
\end{equation*}
$$

with:

$$
\begin{equation*}
\begin{aligned}
\rho_1 &=  \text{min}\left(\frac{(\eta^n)^2}{(\sigma^n)(1+\frac{\alpha}{2})^2\text{TOL}^2},2\right) \\
\rho_2 &=  \text{max}\left(\frac{(\eta^n)^2}{(\sigma^n)(1-\frac{\alpha}{2})^2\text{TOL}^2},0.5\right)
\end{aligned}
\end{equation*}
$$

This choice of $\rho_1$ and $\rho_2$ ensures that the time step size is not increased or decreased by more than a factor of 2 at each time step, in order to avoid too large variations in the time step size, which could lead to instabilities in the simulation. Additionally, the time step size is bounded by a minimum and a maximum value, which are set to $\Delta t_{\text{min}}=10^{-4}$ and $\Delta t_{\text{max}}=0.1$ in this project.

For more details on the specific implementation of the adaptive time-stepping algorithm in this project, please refer to the `adaptive_NS.cc` file, where the code is commented and explained in detail.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Usage

The program can be easily compiled and executed using the `CMake` build system. The following steps are required to compile and run the program:

1. Compile the program using `CMake`: in the main directory of the project, create a build directory and run `CMake` to configure the project:

```bash
cmake -B build -S .
```

Then, after entering the `build` directory, compile the project using `make`, specifying the number of processors `N` to use for the compilation:

```bash
cd build
make -jN
```

2. Move the `nsbench.inp` and `parameter-file.prm` files to the build directory. The `nsbench.inp` file contains the geometry of the test case, while the `parameter-file.prm` file contains the parameters for the simulations. Both files are required for the program to run correctly.
The parameters in the `parameter-file.prm` file can be modified to change the settings of the simulations, such as the Reynolds number, the initial time step size, the initial and final time, the output frequency, the number of refinement levels for the mesh and many others (see the comments in the file for more details).

1. Run the program: to run the program, simply execute the `adaptive_NS` executable in the build directory:

```bash
./adaptive_NS
```

The program will start the simulation and progressively output the results; when running, the code will output the current step of the simulation, the current time and the current time step size. The results of the simulations include the velocity and pressure fields and the vorticity field every `output_interval` time steps, as specified in the `parameter-file.prm` file. The program will also save in the `time_steps.csv` file the time step size evolution during the simulation, which can be later analyzed and plotted using the `time-step_analysis.ipynb` Jupyter notebook.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Results

In this section, the results of some simulations are presented and discussed. The simulations were performed using the adaptive time-stepping algorithm implemented in this project and described in the previous sections <a href="#adaptive-time-stepping-algorithm">here</a>. The objective of these simulations was to analyze the behavior of the adaptive time-stepping algorithm applied to the employed projection method and its impact on the accuracy and the efficiency of the simulations. Finally, the simulations were performed for different Reynolds numbers and different refinement levels of the mesh, in order to study the behavior of the adaptive time-stepping algorithm in different scenarios.

### Test Case

The test case that has been used for the simulations consists of the **flow around a square obstacle in a channel**, which is actually the same presented in the Step-35 tutorial of the deal.II library <a href="#ref1">[1]</a>., which was used as a reference for the development of this project. The geometry of the test case is defined in the `nsbench.inp` file and is shown in the following figure:

<p align="center">
    <img src="images/test_geometry.png" alt="Geometry of the test case" width="700">
</p>

No slip boundary conditions were applied on both the top and bottom walls of the channel and the obstacle, while on the left side of the channel the following inflow boundary condition was imposed:

$$
u = \begin{cases}
    4U_my\frac{(H-y)}{H^2} \\
    0
\end{cases}
$$

with $U_m=1.5$ and $H=4.1$ (the height of the channel). Finally, on the right side of the channel, the condition that the vertical velocity component and the pressure should both be zero was imposed.

As mentioned before, the simulations were performed for **different refinement levels** of the mesh: specifically, the mesh was subjected to two refinement levels, which are denoted as `level1` and `level2`. The `level1` mesh is the coarsest mesh ($960$ active cells), where $2$ uniform refinements (`refine_global` function of deal.II) were applied to the initial mesh, while the `level2` mesh is the finest mesh ($3840$ active cells), where $3$ uniform refinements were applied to the initial mesh. The following figures show the `level1` and `level2` meshes:


<p align="center">
    <img src="images/2_ref_levels.png" alt="Mesh level 1" width="700">
    <p align="center">Level 1 Refinement</p>
</p>

<p align="center">
    <img src="images/3_ref_levels.png" alt="Mesh level 2" width="700">
    <p align="center">Level 2 Refinement</p>
</p>

The simulations were performed for **different Reynolds numbers**, which were set to $Re=50$, $Re=100$ and $Re=200$. The simulations were run for a final time of $T=40$ and the initial time step size was set to $\Delta t_0=0.001$. The output frequency was set to `output_interval=20`, which means that the results were saved every $20$ time steps. The rest of the parameters were set to the values specified in the `parameter-file.prm` file.

### Simulation Results

Here we present the results of the simulation performed for the `level2` mesh for the Reynolds number $Re=100$. The following animation shows the evolution of the velocity magnitude field during the simulation:

<p align="center">
    <img src="images/sol.gif" alt="Velocity magnitude field" width="700">
</p>

The development and extension of a vortex chain behind the obstacle can be observed: the flow is characterized by the formation of vortices and the shedding of vortices downstream of the obstacle. The simulation captures the main features of the flow around the obstacle and the evolution of the vortices in the channel. To better visualize these vortices and the flow structures, the following animation shows the velocity field with a *Line Integral Convolution (LIC)* filter applied:

<p align="center">
    <img src="images/sol_LIC.gif" alt="Velocity field with LIC filter" width="700">
</p>

The LIC filter enhances the visualization of the flow structures and the vortices, making them more visible and highlighting the dynamics of the flow.

The `data` folder contains the results of all the simulations performed for the different Reynolds numbers and refinement levels of the mesh. The results include the velocity and pressure fields and the vorticity field at different time steps of the simulation. The results were saved in the `vtk` format and can be visualized using a software that supports this format, such as `ParaView`. For each simulation, the `time_steps.csv` file containing the time step size evolution during the simulation is also provided.

### Adaptive Time-Stepping Analysis

The `time-step_analysis.ipynb` Jupyter notebook was used to analyze the time step size evolution during the simulations. The notebook reads the different `time_steps.csv` files and plots the time step size evolution over time, as well as over the number of time steps, for each simulation.

We will first provide some plots and then discuss the results of the analysis.

The following figure shows the time step size evolution during two simulations performed with Reynolds number $Re=100$ and the `level1` and `level2` meshes respectively:

<p align="center">
    <img src="images/Ref_Level_Comp.png" alt="Time step size evolution" width="700">
</p>

On the other hand, the following figure shows the time step size evolution during the simulations performed for the different Reynolds numbers and the `level2` mesh, specifically for $Re=50$, $Re=100$ and $Re=200$:

<p align="center">
    <img src="images/Re_Comp.png" alt="Time step size evolution" width="700">
</p>

#### General Observations

It can be noticed that the time step size evolution is different for the different Reynolds numbers and refinement levels of the mesh. However, a general time step behavior can be observed: 

- initially, the time step size decreases, reflecting the initialization of the simulation where transient effects are present and the solution is evolving rapidly; this decrease is not very evident in these plots, but can be observed in the results of the simulations in the `time_steps.csv` files, as well as in the *log* plot of the time step size evolution provided in the notebook. This first short phase corresponds to the need to adapt the time step size to capture the rapid evolution in the flow dynamics;
- as the simulation progresses, the time step size increases, reflecting the stabilization of the flow and its transition to a developing state: the transient effects are reduced and the flow becomes more stable, allowing for larger time step sizes to be used;
- as vortex shedding occurs and the flow structures start to form, the time step size decreases again, reflecting the need to capture the rapid changes in the solution and the flow dynamics; the time step size subsequently oscillates around a certain value with some drops, reflecting the formation and propagation of the vortices and the flow structures;
- finally, the time step size stabilizes into a sort of periodic trend as vortices form and propagate downstream. The small oscillations in the time step size become more regular, mirroring the periodicity of the flow structures and the whole physical flow.
  
Overall, the adaptive time-stepping algorithm successfully adjusts the time step size during the simulation, following the time evolution of the solution and the physical flow dynamics.

#### Effect of the Grid Refinement

The effect of the grid refinement on the time step size evolution can be observed in the first plot: the time step size is generally smaller for the `level2` mesh than for the `level1` mesh, reflecting the need to capture the flow structures and the vortices with a higher resolution when using a finer grid. In the case of the `level2` mesh, the smaller time step size reflects the increased sensitivity of the simulation to gradients and vortex structures, which require a finer mesh to be resolved accurately. On the other hand, the `level1` coarser mesh allows for larger time step sizes; moreover, the trend of the time step size evolution is smoother, more stable and less fluctuating than for the `level2` mesh, indicating the lower resolution of the mesh and the reduced sensitivity of the simulation to the flow structures, which leads to a more stable and less oscillatory behavior of the time step size.

#### Effect of the Reynolds Number

The effect of the Reynolds number on the time step size evolution can be observed in the second plot: the time step size is generally smaller for higher Reynolds numbers. In fact, as $Re$ increases, the flow becomes more turbulent and the vortices form and propagate more rapidly, requiring a smaller time step size to capture the rapid changes in the solution and the flow dynamics. The time step size is larger for lower Reynolds numbers, reflecting the smoother and less turbulent flow, which evolves more slowly and allows for a larger time step size to be resolved accurately enough. The increased oscillatory behavior of the time step size which can be observed for higher Reynolds numbers aligns with the physical flow exhibiting more pronounced vortex shedding, which introduces greater variability in the flow dynamics and the solution over time. Overall, the periodicity of the time step size evolution aligns with the vortex shedding frequency, which is higher for higher Reynolds numbers and lower for lower Reynolds numbers.

### Consistency with the Expected Behavior

All these observations are consistent with the expected behavior of the adaptive time-stepping algorithm, confirming its effectiveness in adjusting the time step size during the simulation to ensure accuracy and efficiency. In fact, the time step size is dynamically modified reflecting the time evolution of the physical flow: it is reduced when the solution is evolving rapidly and the flow structures are forming, and it is increased when the flow is more stable and the solution is changing less, optimizing the computational efficiency of the simulation, without compromising the accuracy of the results. As expected, finer grid refinements and higher Reynolds numbers require smaller time step sizes, due to the increased sensitivity to flow details and the faster dynamics, respectively.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Conclusions

In this project, an **adaptive time-stepping algorithm** for solving the incompressible Navier-Stokes equations with a finite element method using the deal.II library was implemented. The adaptive time-stepping algorithm effectively adjusts the time step size at each time step of the simulation based on a time error indicator, which is used to determine the optimal time step size for the next time step. The time step size is dynamically modified during the simulations, following the time evolution of the solution and the physical flow dynamics, ensuring accuracy and efficiency of the simulations. 

The adaptive time-stepping algorithm was successfully implemented and tested for different Reynolds numbers and refinement levels of the mesh, providing evidence of its effectiveness in controlling the time step size and optimizing the computational efficiency of the simulations. Finally, the results of the simulations were analyzed and discussed, confirming the consistency of the adaptive time-stepping algorithm with the expected behavior and its ability to adapt the time step size to the time evolution of the solution and the physical flow.

In conclusion, it's important to underline the fact that the analysis of the computed time step size during the simulations can help to better identify different phases of the development of the flow: the time step size evolution can provide insights into the flow dynamics and the formation of the flow structures.

Overall, the adaptive time-stepping algorithm can be a valuable tool for controlling the time step size and optimizing the computational efficiency of the simulations, ensuring that the solution is accurate and the computational costs are optimized.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

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
