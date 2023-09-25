# Test that Baby Jubjub represented in Short Weierstrass form is convertible to Montgomery form
# Reference: Conversion criteria in https://en.wikipedia.org/wiki/Montgomery_curve

# Base field
p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
Fr = GF(p)

# Short Weierstrass parameters
SWa = Fr(7296080957279758407415468581752425029516121466805344781232734728849116493472)
SWb = Fr(16213513238399463127589930181672055621146936592900766180517188641980520820846)

# Elliptic curve order must be divisible by 4
ec = EllipticCurve(Fr, [SWa, SWb])
assert(ec.order() % 4 == 0)

# x^3 + SWa * x + SWb must have at least one root in Fr
x = PolynomialRing(Fr, 'x').gen()
f = x^3 + SWa * x + SWb
roots = f.roots()
assert(len(roots) > 0)

# 3 * alpha^2 + SWa must be a quadratic residue in Fr
alpha = Fr(7296080957279758407415468581752425029516121466805344781232734728858602888105)
residue = 3 * alpha^2  + SWa
assert(kronecker(residue, p) == 1)

# Compute the value of the parameter s
s = 1 / residue.sqrt()
print(s)