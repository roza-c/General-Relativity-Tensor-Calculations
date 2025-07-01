# General Relativity Tensor Calculations in Python

This repository contains a Python script that uses the `SymPy` library to perform symbolic tensor calculations for a given spacetime metric within the framework of General Relativity.

The primary goal of this project is to automate the complex and often tedious calculations of key geometric objects, starting from the metric tensor.

## The Spacetime Metric

The script analyzes the geometry of a 3D spacetime defined by the following line element:

$$ds^2 = L^2 [dt^2 - \cosh^2t(d\theta^2 + \sin^2\theta d\phi^2)]$$

where $L$ is a constant. The script first constructs the metric tensor $g_{\mu\nu}$ and its inverse $g^{\mu\nu}$ from this line element.

## Calculations Performed

The script automates the calculation of several key tensors and scalars:

* **Christoffel Symbols:** Computes the Christoffel symbols of the second kind, which describe the effects of parallel transport in a curved spacetime.
    $$\Gamma^\mu_{\nu\lambda} = \frac{1}{2} g^{\mu\sigma} (\partial_\nu g_{\sigma\lambda} + \partial_\lambda g_{\nu\sigma} - \partial_\sigma g_{\nu\lambda})$$

* **Riemann Curvature Tensor:** Calculates the Riemann tensor, which encodes the full curvature of spacetime.
    $$R^\rho_{\sigma\mu\nu} = \partial_\mu \Gamma^\rho_{\sigma\nu} - \partial_\nu \Gamma^\rho_{\sigma\mu} + \Gamma^\lambda_{\sigma\nu}\Gamma^\rho_{\lambda\mu} - \Gamma^\lambda_{\sigma\mu}\Gamma^\rho_{\lambda\nu}$$

* **Ricci Tensor and Scalar:** Finds the Ricci tensor by contracting the Riemann tensor, and then computes the Ricci scalar by tracing the Ricci tensor.
    $$R_{\mu\nu} = R^\lambda_{\mu\nu\lambda}$$
    $$R = g^{\mu\nu}R_{\mu\nu}$$

* **Einstein Field Equations:** Finally, the script verifies that this spacetime is a solution to the Einstein Field Equations with a cosmological constant, $\Lambda$, and solves for the value of $\Lambda$ in terms of $L$.
    $$R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}R + \Lambda g_{\mu\nu} = 0$$

## How to Run the Script

1.  **Prerequisites:** You must have Python and the `SymPy` library installed. If `SymPy` is not installed, you can add it using pip.
    ```bash
    pip install sympy
    ```

2.  **Execute the script:** Run the Python script from your terminal. Be sure to replace `HW2_Q3Q4.py` with the actual name of your script.
    ```bash
    python HW2_Q3Q4.py
    ```

## Results

Upon execution, the script will print the non-zero components of the Christoffel symbols, the Riemann tensor, and the Ricci tensor. It will also display the calculated Ricci scalar and the final value of the cosmological constant $\Lambda$ that satisfies the Einstein Field Equations for this metric.