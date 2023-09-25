from baby_jubjub import BabyJubjubPoint, SWPoint, MontPoint, TwEdPoint

# Verifies different forms of the Baby Jubjub curve are equivalent
# Curve parameters defined in: https://eips.ethereum.org/EIPS/eip-2494
# Conversion reference: https://www-fourier.univ-grenoble-alpes.fr/mphell/doc-v5/conversion_weierstrass_edwards.html
# Note: There is a typo in the conversion reference. The isomorphism f: (x, y) -> (u, v) from Montgomery points 
# to Short Weierstrass points should have a v coordinate of y / B, not x / B
def main():
    # Base field
    Fr = BabyJubjubPoint.Fr

    # Short Weierstrass parameters
    SWa = SWPoint.a
    SWb = SWPoint.b

    # Montgomery parameters
    MontA = MontPoint.A
    MontB = MontPoint.B

    # Twisted Edwards parameters
    TwEdA = TwEdPoint.A
    TwEdd = TwEdPoint.d

    # Check that the curve parameters represent the same curve
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

    # Check that the generator and base points are the same
    assert(SWPoint.generator().to_montgomery() == MontPoint.generator())
    assert(SWPoint.base().to_montgomery() == MontPoint.base())
    assert(MontPoint.generator().to_twisted_edwards() == TwEdPoint.generator())
    assert(MontPoint.base().to_twisted_edwards() == TwEdPoint.base())
    assert(TwEdPoint.generator().to_montgomery() == MontPoint.generator())
    assert(TwEdPoint.base().to_montgomery() == MontPoint.base())
    assert(MontPoint.generator().to_short_weierstrass() == SWPoint.generator())
    assert(MontPoint.base().to_short_weierstrass() == SWPoint.base())

    # Base point is Generator point scaled by 8
    assert(SWPoint.generator().scalar_mul(8) == SWPoint.base())
    assert(MontPoint.generator().scalar_mul(8) == MontPoint.base())
    assert(TwEdPoint.generator().scalar_mul(8) == TwEdPoint.base())

if __name__ == '__main__':
    main()