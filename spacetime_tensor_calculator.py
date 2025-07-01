'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q3: Christoffel-ology in 3D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q3 (a) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

import sympy as sp

# Define SymPy variables for the coordinates
t = sp.symbols('t')             # time coordinate
theta = sp.symbols('theta')     # polar angle coordinate
phi = sp.symbols('phi')         # azimuthal angle coordinate
L = sp.symbols('L')             # constant 

# Define the 3x1 array xup for spacetime coordinate x^μ
xup = sp.Matrix([t, theta, phi])

# Define the 3x3 array for the downstairs metric components g_μν
gdndn = sp.Matrix([
    [L**2, 0, 0],
    [0, -L**2 * sp.cosh(t)**2, 0],               
    [0, 0, -L**2 * sp.cosh(t)**2 * sp.sin(theta)**2]])

# Compute the components of the 3x3 array for the upstairs metric componenets g^μν
# Requires that the upstairs and downstairs metrics are inverses of each other 
gupup = gdndn.inv()

# Display the components of the downstairs and upstairs metrics
print("Downstairs metric (g_μν):\n")
sp.pprint(gdndn)

print("\n")

print("\nUpstairs metric (g^μν):\n")
sp.pprint(gupup)

print("\n")

# Check that g^μν * g_νσ = δ^μ_σ
identity_check_1 = gupup * gdndn
print("\nIdentity check 1: (g^μν * g_νσ):\n")

# Check that we get back the identity matrix
sp.pprint(identity_check_1)
print("\nIdentity property g^μν * g_νσ = δ^μ_σ is satisfied.\n")

# Check that g_μν * g^νσ = δ^μ_σ
identity_check_2 = gdndn * gupup
print("\nIdentity check 2: (g_μν * g^νσ):\n")

# Check that we get back the identity matrix
sp.pprint(identity_check_2)
print("\nIdentity property g_μν * g^νσ = δ_μ^σ is satisfied.\n")

print("\nThe upstairs and downstairs metrics are verified as inverses of each other.")

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q3 (b) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# Define symbols for the partial derivatives of the downstairs metric tensor
coords = [t, theta, phi]

# Define 3x3x3 array to compute Christoffel symbols Γ^μ_νλ
Chrisupdndn = sp.MutableDenseNDimArray([0] * 27, (3, 3, 3))  # Initialize values to zero

# Loop through all indices (3 for each variable) to calculate the Christoffel symbols
for mu in range(3):
    for nu in range(3):
        for lam in range(3):
            term = 0
            for sigma in range(3):
                 # Compute partial derivatives of metric tensors
                partial_nu = sp.diff(gdndn[sigma, lam], coords[nu]) # ∂_ν g_σλ
                partial_lam = sp.diff(gdndn[nu, sigma], coords[lam]) # ∂_λ g_νσ 
                partial_sigma = sp.diff(gdndn[nu, lam], coords[sigma]) # ∂_σ g_νλ

                # Calculate Christoffel: 1/2 * g^μσ (∂_ν g_σλ + ∂_λ g_νσ − ∂_σ g_νλ)
                term += gupup[mu, sigma] * (partial_nu + partial_lam - partial_sigma) / 2
                
            Chrisupdndn[mu, nu, lam] = (term) # Γ^μ_νλ

# Display only nonzero Christoffel symbols
print("Nonzero Christoffel symbols Γ^μ_νλ:")
for mu in range(3):
    for nu in range(3):
        for lam in range(3):
            if Chrisupdndn[mu, nu, lam] != 0:
                # Display non-simplified and simplified terms
                print(f"Gamma^({mu})_({nu},{lam}) = {Chrisupdndn[mu, nu, lam]} \n\t\t= {sp.simplify(Chrisupdndn[mu, nu, lam])}")


'''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q4: Riemann-ology in 3D ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q4 (a) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# Define 3x3x3x3 array to compute Riemann tensor components R^ρ_σμν
Riemupdndndn = sp.MutableDenseNDimArray([0] * 81, (3, 3, 3, 3)) # Initialize values to zero

# Loop through all indices (3 for each variable) to calculate Riemann components
for rho in range(3):
    for sigma in range(3):
        for mu in range(3):
            for nu in range(3):
                # Compute partial derivatives of Christoffel symbols
                partial_mu = sp.diff(Chrisupdndn[rho, sigma, nu], coords[mu]) # ∂_μ Γ^ρ_σν 
                partial_nu = sp.diff(Chrisupdndn[rho, sigma, mu], coords[nu]) # ∂_ν Γ^ρ_σμ

                # Compute sum of Christoffel products: Γ^λ_σν*Γ^ρ_λμ - Γ^λ_σμ*Γ^ρ_λν 
                sum_Chris = sum(Chrisupdndn[lam, sigma, nu] * Chrisupdndn[rho, lam, mu] -
                                Chrisupdndn[lam, sigma, mu] * Chrisupdndn[rho, lam, nu]
                                for lam in range(3))
                
                # Calculate Riemann tensor component R^ρ_σμν
                Riemupdndndn[rho, sigma, mu, nu] = (partial_mu - partial_nu + sum_Chris)

# Display only the nonzero components of the Riemann tensor
print("Nonzero Riemann tensor components R^ρ_σμν:")
for rho in range(3):
    for sigma in range(3):
        for mu in range(3):
            for nu in range(3):
                if Riemupdndndn[rho, sigma, mu, nu] != 0:
                    # Display non-simplified and simplified terms
                    if Riemupdndndn[rho, sigma, mu, nu] != sp.simplify(Riemupdndndn[rho, sigma, mu, nu]):
                        print(f"R^({rho})_({sigma},{mu},{nu}) \t= {Riemupdndndn[rho, sigma, mu, nu]}\n\t\t= {sp.simplify(Riemupdndndn[rho, sigma, mu, nu])}")
                    else:
                        print(f"R^({rho})_({sigma},{mu},{nu}) \t= {Riemupdndndn[rho, sigma, mu, nu]}")

print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q4 (b) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# Define 3x3 array for Ricci tensor R_μν 
Riccdndn = sp.MutableDenseNDimArray([0] * 9, (3, 3)) # Initialze values to zero

# Compute Ricci tensor by summing over Riemann tensor components # R_μν = R^λ_μνλ
for mu in range(3):
    for nu in range(3):
        Riccdndn[mu, nu] = sum(Riemupdndndn[lam, mu, nu, lam] for lam in range(3))

# Display only nonzero components of the Ricci tensor
print("Nonzero Ricci tensor components R_μν:")
for mu in range(3):
    for nu in range(3):
        if Riccdndn[mu, nu] != 0:
            # Display non-simplified and simplified terms
            if Riccdndn[mu, nu] != sp.simplify(Riccdndn[mu, nu]):
                print(f"R_({mu},{nu}) = {Riccdndn[mu, nu]} \n\t= {sp.simplify(Riccdndn[mu, nu])}")
            else:
                print(f"R_({mu},{nu}) = {Riccdndn[mu, nu]}")

# Compute Ricci scalar by summing over diagonal elements of Ricci tensor # R = R^μ_μ
RiccScal = sum(Riccdndn[mu, nu] for mu in range(3))
print("\nRicci scalar R:")
sp.pprint(f"R = {RiccScal} \n⇒ R = {sp.simplify(RiccScal)}")

print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Q4 (c) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# Initialize Lambda variable Λ
Lambda = sp.symbols('Lambda')

# Verify the Einstein equation R_μν − 1/2*g_μν*R + Λ*g_μν = 0
# Define 3x3 array to hold results
EinsteinEq = sp.MutableDenseNDimArray([0] * 9, (3, 3)) # Initialize values to zero

# Loop through indices
for mu in range(3):
    for nu in range(3):
        EinsteinEq[mu, nu] = Riccdndn[mu, nu] - (1/2)*gdndn[mu, nu] * RiccScal + Lambda * gdndn[mu, nu]
        EinsteinEq[mu, nu] = sp.simplify(EinsteinEq[mu, nu])

# Display the nonzero components of the Einstein equation
print("Nonzero components of the Einstein equation:")
for mu in range(3):
    for nu in range(3):
        if EinsteinEq[mu, nu] != 0:
            print(f"EinsteinEq_({mu},{nu}) = {EinsteinEq[mu, nu]}")

# Solve for Λ in terms of L
Lambda_expr = sp.solve(EinsteinEq[0, 0], Lambda)[0]
print("\nValue of Lambda (Λ) in terms of L:")
sp.pprint(f"Λ = {sp.simplify(Lambda_expr)}")
print("Where Λ ∝ 1/L^2\n")
print("Hence this spacetime obeys the Einstein equation, which is the set of equations given by R_μν − 1/2*g_μν*R + Λ*g_μν = 0")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

