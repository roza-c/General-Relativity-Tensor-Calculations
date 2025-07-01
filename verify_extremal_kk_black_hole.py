import sympy as sp

# Define SymPy variables for the coordinates and constants
t, r, theta, phi = sp.symbols('t r theta phi', real=True, positive=True)
q = sp.Symbol('q', real=True, positive=True)
H = 1 + q/r
chi_expr = -sp.sqrt(3)/2 * sp.log(H) 
exp_fac = sp.exp(-sp.sqrt(3)*chi_expr)  

# Define 4x4 array for downstairs metric (gd) components 
# Indices: (t=0, r=1, theta=2, phi=3)
gd = sp.Matrix([
    [1/sp.sqrt(H), 0, 0, 0],
    [0, -sp.sqrt(H), 0, 0],               
    [0, 0, -r**2*sp.sqrt(H), 0],
    [0, 0, 0, -r**2*sp.sqrt(H)*(sp.sin(theta))**2]])

print("The 4D Kaluza-Klein Black Hole (downstairs) Metric g_{μν}")
sp.pprint(gd)
print()

# Compute upstairs metric as inverse of downstairs metric
gup = gd.inv()

print("The Inverse (upstairs) Metric g^{μν}")
sp.pprint(gup)
print()

# Define coordinate symbols for indexing
coords = [t, r, theta, phi]

# Initialize 4x4x4 array to compute Christoffel symbols
Gamma = sp.MutableDenseNDimArray.zeros(4,4,4)

# Loop through all indices to calculate Christoffels
for mu in range(4):
    for nu in range(4):
        for lam in range(4):
            # sum over sigma
            term = 0
            for sigma in range(4):
                # Compute partial derivatives of metric tensors
                partial_nu = sp.diff(gd[sigma, lam], coords[nu])
                partial_lam = sp.diff(gd[nu, sigma], coords[lam])
                partial_sigma = sp.diff(gd[nu, lam], coords[sigma])
                term += gup[mu, sigma] * (partial_nu + partial_lam - partial_sigma) / 2
            Gamma[mu, nu, lam] = sp.simplify(term)

print("Nonzero Christoffel Symbols Γ^μ_(νλ):")
for mu in range(4):
    for nu in range(4):
        for lam in range(4):
            if (Gamma[mu, nu, lam]) != 0:
                print(f"Gamma^({mu})_({nu},{lam}) = {Gamma[mu, nu, lam]}")
print()

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

                # Compute sum of Christoffel products over lambda
                sum = 0
                for lam in range(4):
                    sum += Gamma[rho, mu, lam]*Gamma[lam, nu, sigma] - Gamma[rho, nu, lam]*Gamma[lam, mu, sigma]
                Riem[rho, sigma, mu, nu] = sp.simplify(term + sum)

print("Nonzero Riemann tensor components R^{ρ}_{σμν}:")
for rho in range(4):
    for sigma in range(4):
        for mu in range(4):
            for nu in range(4):
                if Riem[rho, sigma, mu, nu] != 0:
                    print(f"R^({rho})_({sigma},{mu},{nu}) \t= {sp.simplify(Riem[rho, sigma, mu, nu])}")

# Define 4x4 array for Ricci tensor
Ric = sp.MutableDenseNDimArray.zeros(4,4)

# Compute Ricci by summing over Riemann tensor components
for mu in range(4):
    for nu in range(4):
        val = 0
        for rho in range(4):
            val += Riem[rho, mu, rho, nu]
        Ric[mu, nu] = sp.simplify(val)

print("\nNonzero Ricci tensor components R_{μν}:")
for mu in range(4):
    for nu in range(4):
        if Ric[mu, nu] != 0:
            print(f"R_({mu},{nu}) = {sp.simplify(Ric[mu, nu])}")

# Compute Ricci scalar by summing over diagonal elements of Ricci tensor
RicScal = 0
for mu in range(4):
    for nu in range(4):
        RicScal += gup[mu, nu]*Ric[mu, nu]
RicScal = sp.simplify(RicScal)

print("\nRicci Scalar R:")
sp.pprint(RicScal)
print()

# Compite partial derivatives of chi wrt t,r,theta,phi
# [dchi/dt, dchi/dr, dchi/dθ, dchi/dφ]
dchi = [sp.diff(chi_expr, x) for x in coords]  

# Compute (\partial\chi)^2 = g^{\mu\nu} (\partial\mu\chi)(\partial\nu\chi)
dchi_sq = 0
for mu in range(4):
    for nu in range(4):
        dchi_sq += gup[mu, nu]*dchi[mu]*dchi[nu]
dchi_sq_simpl = sp.simplify(dchi_sq)

# Compute T_{\mu\nu}[\chi] = \partial\mu\chi \partial\nu\chi - 1/2 g_{\mu\nu} (\partial\chi)^2
Tchi = sp.MutableDenseNDimArray.zeros(4,4)
for mu in range(4):
    for nu in range(4):
        val = dchi[mu]*dchi[nu] - sp.Rational(1,2)*gd[mu, nu]*dchi_sq_simpl
        Tchi[mu, nu] = sp.simplify(val)

# Define A\mu: A0(t,r,theta,phi) = 1 - 1/H(r); others=0
A = [ (1 - 1/H), 0, 0, 0 ]

# Compute F_{\mu\nu}
Fdn = sp.MutableDenseNDimArray.zeros(4,4)
for mu in range(4):
    for nu in range(4):
        Fdn[mu, nu] = sp.diff(A[nu], coords[mu]) - sp.diff(A[mu], coords[nu])
Fdn_simpl = [[sp.simplify(Fdn[mu,nu]) for nu in range(4)] for mu in range(4)]

# Raise indices: F^\mu\nu = g^{\mu\alpha} g^{\nu\beta} F_{\alpha\beta}
Fup = sp.MutableDenseNDimArray.zeros(4,4)
for mu in range(4):
    for nu in range(4):
        val = 0
        for alpha in range(4):
            for beta in range(4):
                val += gup[mu, alpha]*gup[nu, beta]*Fdn[alpha, beta]
        Fup[mu, nu] = sp.simplify(val)

# Contract F^2 = F^{\mu\nu} F_{\mu\nu}
F_squared = 0
for mu in range(4):
    for nu in range(4):
        F_squared += Fup[mu, nu]*Fdn[mu, nu]
F_squared_simpl = sp.simplify(F_squared)

# First compute F\mu^\lambda = g^{\lambda\sigma} F_{\mu\sigma}
Fmuldn_up = sp.MutableDenseNDimArray.zeros(4,4)
for mu in range(4):
    for lam in range(4):
        val = 0
        for sigma in range(4):
            val += gup[lam, sigma]*Fdn[mu, sigma]
        Fmuldn_up[mu, lam] = sp.simplify(val)

# Then build T_{\mu\nu}[A]
TA = sp.MutableDenseNDimArray.zeros(4,4)
for mu in range(4):
    for nu in range(4):
        # (F\mu^\lambda F\nu\lambda) 
        Fmu_lam_FnuLam = 0
        for lam in range(4):
            Fmu_lam_FnuLam += Fmuldn_up[mu, lam]*Fdn[nu, lam]
        val = -exp_fac*( Fmu_lam_FnuLam - sp.Rational(1,4)*gd[mu, nu]*F_squared_simpl )
        TA[mu, nu] = sp.simplify(val)

# Check if LHS - RHS = 0
Check = sp.MutableDenseNDimArray.zeros(4,4)
for mu in range(4):
    for nu in range(4):
        # LHS = R_{\mu\nu} - 1/2 g_{\mu\nu}R
        lhs = Ric[mu, nu] - sp.Rational(1,2)*gd[mu, nu]*RicScal
        # RHS = -1/2 T_{\mu\nu}
        rhs = -sp.Rational(1,2)*(Tchi[mu, nu] + TA[mu, nu])
        eqn = lhs + rhs
        Check[mu, nu] = sp.simplify(eqn)

print("Einstein Equation Check: R_{μν} - 1/2*g_{μν}*R + 1/2*(T_{μν}[χ] + T_{μν}[A])")
nonzero_flag = False
for mu in range(4):
    for nu in range(4):
        val = sp.simplify(Check[mu, nu])
        if val != 0:
            nonzero_flag = True
            print(f"Check({mu},{nu}) = {val}")
if not nonzero_flag:
    print("All components vanish, thus the Einstein EOM is satisfied!")
