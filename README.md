# Advanced 3D Graphics & Physics Engine Portfolio
**Ruyi (Evelyn) Bao | M.S. in Computer Science @ USC**

This repository features a suite of high-performance C++ modules focusing on low-level computer graphics, numerical physics, and skeletal animation. All solvers and rendering pipelines were implemented from scratch.

---

## 🎥 Project Demos
### Animation & Physics (520 Series)
| Inverse Kinematics (IK) | MoCap Interpolation | Soft Body Physics |
| :---: | :---: | :---: |
| ![IK Demo](./Demo/520-3.gif) | ![MoCap Demo](./Demo/520-2-bq.gif) | ![Jello Demo](./Demo/520-1.gif) |
| *Jacobian Solver & Regularization* | *Quaternion SLERP & Animation*(Note: High-fidelity GIF, may load slowly) | *Mass-Spring System (RK4)* |

### Rendering & Geometry (420 Series)
| Roller Coaster Simulation | Recursive Raytracing |
| :---: | :---: |
| <video src="./Demo/420-2.mp4" width="100%" muted autoplay loop></video> | ![Raytrace Result](./Rendering-Geometry/420-3-Raytrace/toy.jpg) |
| *Catmull-Rom Splines (Click for MP4)* | *Phong Reflection & Hard Shadows* |
---

## 🛠️ Core Technical Highlights

### 1. Animation & Simulation (CSCI 520)
* **Jacobian IK Solver**: Implemented a robust IK engine using **Jacobian matrices** and **ADOL-C** for automatic differentiation. Used **Tikhonov Regularization** to maintain stability near singularities.
* **Soft Body Dynamics**: A 3D "Jello Cube" simulation based on a **Mass-Spring system**. To ensure numerical stability, I implemented a **4th-order Runge-Kutta (RK4)** integrator.
* **MoCap Processing**: Developed a parser for **ASF/AMC** data. Features **Quaternion SLERP** and **Cubic Bezier Splines** to ensure smooth, $C^1$ continuous motion.



### 2. Rendering & Geometry (CSCI 420)
* **Recursive Ray Tracer**: Hand-coded intersection logic for spheres and triangles. Supports recursive reflections and **Phong lighting**.
* **Spline-Based Simulation**: Generated 3D rollercoaster tracks using **Catmull-Rom Splines**. Integrated **Frenet-Serret frames** to derive consistent camera orientation along the path.
* **GLSL Terrain**: Optimized heightmap rendering utilizing modern OpenGL **VBO/VAO** and custom **Vertex/Fragment shaders**.

---

## 📐 Mathematical Foundations

* **Numerical Integration (RK4)** for stability:
$$\mathbf{y}_{n+1} = \mathbf{y}_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

* **Damped Least Squares (IK Stability)**:
$$\Delta \theta = J^T (JJ^T + \lambda^2 I)^{-1} \vec{e}$$

---

## 💻 Tech Stack
* **Languages**: C++ (11/14), GLSL
* **Libraries**: Eigen, ADOL-C, OpenGL, GLUT
* **Build System**: Custom Makefiles

---

## 📜 Academic Context
Developed during the **CSCI 420 (Computer Graphics)** and **CSCI 520 (Computer Animation)** courses at the **University of Southern California**.