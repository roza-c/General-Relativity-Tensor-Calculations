import sympy as sp

# Define SymPy variables, coordinates and functions
t, r, theta, phi = sp.symbols('t r theta phi')
mu, a, Lambda = sp.symbols('mu a Lambda')
A = 1 - 2*mu/r - (Lambda*r**2)/3
alpha = sp.log(A)/2
beta  = -sp.log(A)/2

# Define coordinate symbols for indexing
coords = [t, r, theta, phi]

# Define downstairs metric 
gd = sp.Matrix([
    [sp.exp(2*alpha), 0, 0, 0],
    [0, -sp.exp(2*beta), 0, 0],               
    [0, 0, -r**2, 0],
    [0, 0, 0, -r**2*(sp.sin(theta))**2]])

print("\nDownstairs metric g_{μ,ν}:\n")
sp.pprint(gd)

# Compute upstairs metric as inverse of downstairs metric
gup = gd.inv()

print("\nUpstairs metric g^{μ,ν}:\n")
sp.pprint(gup)

# Initialize 4x4x4 array to compute Christoffel symbols
Gamma = sp.MutableDenseNDimArray.zeros(4,4,4)

# Loop through indices to calculate Christoffels
for mu in range(4):
    for nu in range(4):
        for lam in range(4):
            term = 0 # Sum over sigma
            for sigma in range(4):
                # Compute partial derivatives of metric tensors
                partial_nu = sp.diff(gd[sigma, lam], coords[nu])
                partial_lam = sp.diff(gd[nu, sigma], coords[lam])
                partial_sigma = sp.diff(gd[nu, lam], coords[sigma])
                term += gup[mu, sigma]/2 * (partial_nu + partial_lam - partial_sigma)
            Gamma[mu, nu, lam] = sp.simplify(term)

print("\nNonzero Christoffel Symbols Γ^μ_{νλ}:")
for mu in range(4):
    for nu in range(4):
        for lam in range(4):
            if (Gamma[mu, nu, lam]) != 0:
                print(f"Gamma^({mu})_({nu},{lam}) = {Gamma[mu, nu, lam]}")

# Initialize 4x4x4x4 array to compute Riemann tensor components
Riem = sp.MutableDenseNDimArray.zeros(4,4,4,4)
for rho in range(4):
    for sigma in range(4):
        for mu in range(4):
            for nu in range(4):
                # Compute partial derivatives of Christoffels
                partial_mu = sp.diff(Gamma[rho, sigma, nu], coords[mu])
                partial_nu = sp.diff(Gamma[rho, sigma, mu], coords[nu]) 
                term = partial_mu - partial_nu
                sum = 0 # Sum over lambda
                for lam in range(4): 
                    sum += Gamma[lam, sigma, nu]*Gamma[rho, lam, mu] - Gamma[lam, sigma, mu]*Gamma[rho, lam, nu]
                Riem[rho, sigma, mu, nu] = sp.simplify(term + sum)

print("\nNonzero Riemann tensor components R^{ρ}_{σμν}:")
for rho in range(4):
    for sigma in range(4):
        for mu in range(4):
            for nu in range(4):
                if Riem[rho, sigma, mu, nu] != 0:
                    print(f"R^({rho})_({sigma},{mu},{nu}) \t= {(Riem[rho, sigma, mu, nu])}")

# Define 4x4 array for Ricci tensor
Ric = sp.MutableDenseNDimArray.zeros(4,4)

# Compute Ricci by summing over Riemann tensor components
for mu in range(4):
    for nu in range(4):
        val = 0 # Sum over rho
        for rho in range(4):
            val += Riem[rho, mu, nu, rho]
        Ric[mu, nu] = sp.simplify(val)

print("\nNonzero Ricci tensor components R_{μν}:")
for mu in range(4):
    for nu in range(4):
        if Ric[mu, nu] != 0:
            print(f"R_({mu},{nu}) = {(Ric[mu, nu])}")

# Compute Ricci scalar by summing over diagonal elements of Ricci tensor
RicScal = 0
for mu in range(4):
    for nu in range(4):
        RicScal += gup[mu, nu]*Ric[mu, nu]
RicScal = sp.simplify(RicScal)

print("\nRicci Scalar R:")
sp.pprint(RicScal)

# Compute f(R), f'(R)
f_R = RicScal + a*RicScal**2 - 2*Lambda
fprime_R = 1 + 2*a*RicScal

# Compute Del f'(R)
partial_fprime = [0,0,0,0]
for nu in range(4):
    partial_fprime[nu] = sp.simplify(sp.diff(fprime_R, coords[nu]))

# Define 4x4 array for the second covariant derivative of f'(R)
# Second covariant derivative of a scalar field is just the partial derivative
del_mu_del_nu_fprime = sp.MutableDenseNDimArray.zeros(4,4)
for mu in range(4):
    for nu in range(4):
        del_mu_del_nu_fprime[mu, nu] = sp.simplify(sp.diff(partial_fprime[nu], coords[mu]))

# Compute box(f'(R))
box_fprime = 0
for mu in range(4):
    for nu in range(4):
        box_fprime += gup[mu, nu]*del_mu_del_nu_fprime[mu, nu]
box_fprime = sp.simplify(box_fprime)

# Define 4x4 array to compute the EOM for f(R) gravity 
LHS = sp.MutableDenseNDimArray.zeros(4,4)

# Compute the LHS
for mu in range(4):
    for nu in range(4):
        term = (fprime_R*Ric[mu, nu] - (gd[mu, nu]*f_R)/2 + gd[mu, nu]*box_fprime - del_mu_del_nu_fprime[mu, nu])
        LHS[mu, nu] = sp.simplify(term)

# Check if any components are zero
print("\nCheck EOM LHS_{μν} = 0")
any_nonzero = False
for mu in range(4):
    for nu in range(4):
        lhs = sp.simplify(LHS[mu, nu])
        if lhs != 0: # Print out the components of LHS that are non-zero
            any_nonzero = True
            print(f"LHS[{mu},{nu}] = {lhs}")

if not any_nonzero:
    sp.pprint(LHS)
    print(f"\nAll components of the EOM are zero!\nHence our cosmological black holes solve the modified gravity field equations.\n")
