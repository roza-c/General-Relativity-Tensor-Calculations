import sympy as sp

# Define SymPy variables, coordinates and functions
t, r, theta, phi = sp.symbols('t r theta phi')
mu, a, Lambda, q, G_N = sp.symbols('mu a Lambda q G_N')
A = 1 - (2*mu)/r + (q**2 * (1 + 8*a*Lambda)) / r**2 - (Lambda * r**2) / 3
alpha = sp.simplify(sp.Rational(1,2)*sp.log(A))
beta = sp.simplify(-sp.Rational(1,2)*sp.log(A))

# Define coordinate symbols for indexing
coords = [t, r, theta, phi]

# Define downstairs metric 
gd = sp.simplify(sp.Matrix([
    [sp.exp(2*alpha), 0, 0, 0],
    [0, -sp.exp(2*beta), 0, 0],               
    [0, 0, -r**2, 0],
    [0, 0, 0, -r**2*(sp.sin(theta))**2]]))

print("\nDownstairs metric g_{μ,ν}:\n")
sp.pprint(gd)

# Compute upstairs metric as inverse of downstairs metric
gup = gd.inv()

print("\nUpstairs metric g^{μ,ν}:\n")
sp.pprint(gup)

# Initialize 4x4x4 array to compute Christoffel symbols
Gamma = sp.MutableDenseNDimArray.zeros(4,4,4)

# Loop through indices to calculate Christoffels
for mu_i in range(4):
    for nu in range(4):
        for lam in range(4):
            term = 0 # Sum over sigma
            for sigma in range(4):
                # Compute partial derivatives of metric tensors
                partial_nu = sp.diff(gd[sigma, lam], coords[nu])
                partial_lam = sp.diff(gd[nu, sigma], coords[lam])
                partial_sigma = sp.diff(gd[nu, lam], coords[sigma])
                term += sp.Rational(1,2)*gup[mu_i, sigma] * (partial_nu + partial_lam - partial_sigma)
            Gamma[mu_i, nu, lam] = sp.simplify(term)
            
print("\nNonzero Christoffel Symbols Γ^μ_{νλ}:")
for mu_i in range(3):
    for nu in range(3):
        for lam in range(3):
            if Gamma[mu_i, nu, lam] != 0:
                # Display non-simplified and simplified terms
                print(f"Gamma^({mu_i})_({nu},{lam}) = {(Gamma[mu_i, nu, lam])}")

# Initialize 4x4x4x4 array to compute Riemann tensor components
Riem = sp.MutableDenseNDimArray.zeros(4,4,4,4)
for rho in range(4):
    for sigma in range(4):
        for mu_i in range(4):
            for nu in range(4):
                # Compute partial derivatives of Christoffels
                partial_mu = sp.diff(Gamma[rho, sigma, nu], coords[mu_i])
                partial_nu = sp.diff(Gamma[rho, sigma, mu_i], coords[nu]) 
                term = partial_mu - partial_nu
                sum = 0 # Sum over lambda
                for lam in range(4): 
                    sum += Gamma[lam, sigma, nu]*Gamma[rho, lam, mu_i] - Gamma[lam, sigma, mu_i]*Gamma[rho, lam, nu]
                Riem[rho, sigma, mu_i, nu] = sp.simplify(term + sum)

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
for mu_i in range(4):
    for nu in range(4):
        val = 0
        for rho in range(4):
            val += Riem[rho, mu_i, nu, rho]  
        Ric[mu_i, nu] = sp.simplify(val)

print("\nNonzero Ricci tensor components R_{μν}:")
for i in range(4):
    for j in range(4):
        if Ric[i, j] != 0:
            print(f"R_({i},{j}) = {(Ric[i, j])}")

# Compute Ricci scalar by summing over diagonal elements of Ricci tensor
RicScal = 0
for mu_i in range(4):
    for nu in range(4):
        RicScal += gup[mu_i, nu]*Ric[mu_i, nu]
RicScal = sp.simplify(RicScal)

print("\nRicci Scalar R:")
sp.pprint(RicScal)

# Compute LHS from Q4
f_R = RicScal + a*RicScal**2 - 2*Lambda
fprime_R = 1 + 2*a*RicScal

# Compute Del f'(R)
partial_fprime = [0,0,0,0]
for nu in range(4):
    partial_fprime[nu] = sp.simplify(sp.diff(fprime_R, coords[nu]))

# Define 4x4 array for the second covariant derivative of f'(R)
del_mu_del_nu_fprime = sp.MutableDenseNDimArray.zeros(4,4)
for mu_i in range(4):
    for nu in range(4):
        del_mu_del_nu_fprime[mu_i, nu] = sp.simplify(sp.diff(partial_fprime[nu], coords[mu_i]))

# Compute box(f'(R))
box_fprime = 0
for mu_i in range(4):
    for nu in range(4):
        box_fprime += gup[mu_i, nu]*del_mu_del_nu_fprime[mu_i, nu]
box_fprime = sp.simplify(box_fprime)

# Define 4-potential A_mu(R)
A_dn = sp.zeros(4, 1)
A_dn[0] = q**2/r  

# Define 4x4 array to compute the field strength F_{mu,nu} (downstairs)
F_dndn = sp.MutableDenseNDimArray.zeros(4,4)
for mu_i in range(4):
    for nu in range(4): # = partial_mu A_nu - partial_nu A_mu
        term1 = sp.diff(A_dn[nu], coords[mu_i])
        term2 = sp.diff(A_dn[mu_i], coords[nu])
        F_dndn[mu_i, nu] = sp.simplify(term1 - term2)

print("\nF_dndn is:")
sp.pprint(F_dndn)

# Define 4x4 array to compute field strength F_{nu}^{alpha}
# F_{nu}^{alpha} = F_{nu,beta}g^{beta,alpha}  
F_dn_up = sp.MutableDenseNDimArray.zeros(4,4)
for nu in range(4):
    for alpha in range(4):
        val = 0
        for beta in range(4):
            val += F_dndn[nu, beta] * gup[beta, alpha]
        F_dn_up[alpha, nu] = sp.simplify(val)

print("\nF_dnup is:")
sp.pprint(F_dn_up)

# Define 4x4 array to compute the field strength F^{mu,nu} (upstairs)
# F^{alpha,beta} = g^{alpha,mu} g^{beta,nu} F_{mu,nu}
F_upup = sp.MutableDenseNDimArray.zeros(4,4)
for alpha in range(4):
    for beta in range(4):
        val = 0
        for mu_i in range(4):
            for nu in range(4):
                val += gup[alpha, mu_i]* gup[beta, nu] * F_dndn[mu_i, nu]
        F_upup[alpha, beta] = sp.simplify(val)

# Compute F^{alpha,beta}F_{alpha,beta}
F2 = 0
for alpha in range(4):
    for beta in range(4):
        F2 += F_upup[alpha,beta]*F_dndn[alpha,beta]
F2 = sp.simplify(F2)

print("\nF^2 is:")
sp.pprint(F2)

# Compute EM energy-momentum tensor T_{mu,nu}
# T_{mu,nu} = -(F_{mu alpha}F_{nu}^{alpha}+1/4 g_{mu,nu} F^2)
T_munu = sp.MutableDenseNDimArray.zeros(4,4)
for mu_i in range(4):
    for nu in range(4):
        term1 = 0
        for alpha in range(4):
            term1 += F_dndn[mu_i, alpha]*F_dn_up[nu, alpha]
        term2 = sp.Rational(1,4) * gd[mu_i, nu] * F2
        T_munu[mu_i, nu] = sp.simplify(-term1 + term2)

# Compute and print the LHS
LHS = sp.MutableDenseNDimArray.zeros(4,4)
for mu_i in range(4):
    for nu in range(4):
        term = (fprime_R*Ric[mu_i, nu] - (gd[mu_i, nu]*f_R)/2 + gd[mu_i, nu]*box_fprime - del_mu_del_nu_fprime[mu_i, nu])
        LHS[mu_i, nu] = sp.simplify(term)

print("\n")
for i in range(4):
    for j in range(4):
        print(f"LHS[{i},{j}] =", LHS[i,j].simplify())

# Compute and print the RHS
RHS = sp.MutableDenseNDimArray.zeros(4,4)
for mu_i in range(4):
    for nu in range(4):
        term = -8*sp.pi*G_N*T_munu[mu_i, nu]
        RHS[mu_i, nu] = sp.simplify(term)

print("\n")
for i in range(4):
    for j in range(4):
        print(f"RHS[{i},{j}] =", RHS[i,j].simplify())

# Compute the difference and check for nonzero components
Diff = sp.MutableDenseNDimArray.zeros(4,4)
for mu_i in range(4):
    for nu in range(4):
        Diff[mu_i, nu] = sp.simplify(LHS[mu_i, nu] - RHS[mu_i, nu])

print("\n")
all_zero = True
for mu_i in range(4):
    for nu in range(4):
        if Diff[mu_i, nu] != 0:
            print(f"Diff[{mu_i},{nu}] = {Diff[mu_i,nu]}")
            all_zero = False

if all_zero:
    print("\nAll components vanish!")
else:
    print("\nSome components are nonzero.")
