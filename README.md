# General Relativity Tensor Calculations in Python

This repository contains a collection of Python scripts that use the `SymPy` library to perform symbolic tensor calculations in General Relativity and modified theories of gravity. Each script tackles a specific problem, from calculating the geometric properties of a spacetime to verifying solutions of complex field equations.

## How to Run the Scripts

1.  **Prerequisites:** You must have Python and the `SymPy` library installed. If `SymPy` is not installed, you can add it using `pip`.
    ```bash
    pip install sympy
    ```

2.  **Execute a script:** Run the desired Python script from your terminal. Replace `<script_name.py>` with the name of the file you wish to run.
    ```bash
    python <script_name.py>
    ```

---

## Projects

### Project 1: Curvature of a 3D de Sitter-like Spacetime

* **Script:** `spacetime_tensor_calculator.py` (or similar)

This script analyzes the geometry of a 3D spacetime to calculate its Christoffel symbols, Riemann curvature tensor, Ricci tensor, and Ricci scalar. It verifies that the metric is a solution to the Einstein Field Equations with a cosmological constant.

### Project 2: Extremal Kaluza-Klein Black Hole Solution

* **Script:** `verify_extremal_kk_black_hole.py`

This script verifies that the extremal Kaluza-Klein (KK) black hole metric is a valid solution to the coupled Einstein-Maxwell-Scalar field equations. It computes the stress-energy tensors for the scalar and gauge fields and confirms that the Einstein Field Equations are satisfied.

### Project 3: Schwarzschild-de Sitter Black Hole in f(R) Gravity (Q4)

* **Script:** `verify_fr_schw_desitter.py`

This script demonstrates that a Schwarzschild-de Sitter black hole is a vacuum solution to the field equations of a specific $f(R)$ gravity theory.

* **Theory:** The $f(R)$ model is given by:
    $$f(R) = -2\Lambda + R + aR^2$$

* **Metric:** The Schwarzschild-de Sitter metric is:
    $$ds^2 = \left(1 - \frac{2\mu}{r} - \frac{\Lambda r^2}{3}\right)dt^2 - \left(1 - \frac{2\mu}{r} - \frac{\Lambda r^2}{3}\right)^{-1}dr^2 - r^2 d\Omega^2$$

The script programmatically verifies that this metric and theory satisfy the $f(R)$ equations of motion in the Jordan frame for a vacuum ($T_{\mu\nu}=0$):

$$f'(R)R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}f(R) + [g_{\mu\nu}\Box - \nabla_\mu\nabla_\nu]f'(R) = 0$$

### Project 4: Reissner-Nordström Black Hole in f(R) Gravity (Q5)

* **Script:** `verify_fr_reissner_nordstrom.py`

This script shows that a charged, cosmological black hole is a solution to the same $f(R)$ gravity theory as in Project 3, but now with an electromagnetic field present.

* **Metric and Gauge Field:** The Reissner-Nordström-de Sitter solution is altered in this theory to:
    $$e^{2\alpha(r)} = e^{-2\beta(r)} = 1 - \frac{2\mu}{r} + \frac{q^2(1 + 8a\Lambda)}{r^2} - \frac{\Lambda r^2}{3}$$
    $$A_t = \frac{q^2}{r}$$

The script verifies the full $f(R)$ equations of motion with the electromagnetic stress-energy tensor $T_{\mu\nu}$ as a source:

$$f'(R)R_{\mu\nu} - \frac{1}{2}g_{\mu\nu}f(R) + [g_{\mu\nu}\Box - \nabla_\mu\nabla_\nu]f'(R) = -8\pi G_N T_{\mu\nu}$$

It does this by computing all terms on the left-hand side, computing the electromagnetic stress-energy tensor for the given gauge field $A_t$, and showing that the two sides of the equation are equal.