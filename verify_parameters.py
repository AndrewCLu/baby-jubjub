from algebra import GF
from baby_jubjub import SWPoint

# Verifies different forms of the Baby Jubjub curve are equivalent
# Reference: https://www-fourier.univ-grenoble-alpes.fr/mphell/doc-v5/conversion_weierstrass_edwards.html
def main():
    # Base field
    p=21888242871839275222246405745257275088548364400416034343698204186575808495617
    Fr=GF(p)

    # Short Weierstrass parameters
    SWa=Fr(7296080957279758407415468581752425029516121466805344781232734728849116493472)
    SWb=Fr(16213513238399463127589930181672055621146936592900766180517188641980520820846)

    # Twisted Edwards parameters
    TwEdA = Fr(168700)
    TwEdd = Fr(168696)

    # Montgomery parameters
    MontA = Fr(168698)
    MontB = Fr(1)

    # Convert from Montgomery to Edwards
    assert((MontA + Fr(2)) / MontB == TwEdA)
    assert((MontA - Fr(2)) / MontB == TwEdd)
    
    # Convert from Edwards to Montgomery
    assert(Fr(2) * (TwEdA + TwEdd) / (TwEdA - TwEdd) == MontA)
    assert(Fr(4) / (TwEdA - TwEdd) == MontB)

    # Convert from Montgomery to Weierstrass
    assert((Fr(1) / (MontB * MontB)) * (Fr(1) - (MontA * MontA)/Fr(3)) == SWa)
    assert((MontA / (Fr(3) * MontB * MontB * MontB)) * (Fr(2) * MontA * MontA / Fr(9) - Fr(1)) == SWb)

    # Convert from Weierstrass to Montgomery
    alpha = MontA / Fr(3)
    beta = Fr(1) / MontB
    assert(Fr(3) * alpha * alpha + SWa == beta * beta)
    pt = SWPoint(alpha, 0)
    assert(pt + pt == SWPoint(None, None))

if __name__ == '__main__':
    main()