from baby_jubjub import SWPoint, MontPoint, TwEdPoint

# Verifies different forms of the Baby Jubjub curve are equivalent
# Reference: https://www-fourier.univ-grenoble-alpes.fr/mphell/doc-v5/conversion_weierstrass_edwards.html
def main():
    # Test some operations using Short Weierstrass Generator point
    scalar = 1023
    random_multiple = 1234567890
    num_points_to_test = 100
    for i in range(0, SWPoint.order, SWPoint.order // num_points_to_test):
        print(f"Verifying conversions for {i}G")
        # Generate a point on the curve
        SwG = SWPoint.generator().scalar_mul(i)
        assert(SwG.is_on_curve())

        # Convert from Short Weierstrass to Montgomery
        MontG = SwG.to_montgomery()
        assert(MontG.is_on_curve())

        # Convert from Montgomery to Twisted Edwards
        TwEdG = MontG.to_twisted_edwards()
        assert(TwEdG.is_on_curve())

        # Convert from Twisted Edwards to Montgomery
        MontGPrime = TwEdG.to_montgomery()
        assert(MontGPrime.is_on_curve())
        assert(MontG == MontGPrime)

        # Convert from Montgomery to Short Weierstrass
        SwGPrime = MontGPrime.to_short_weierstrass()
        assert(SwGPrime.is_on_curve())
        assert(SwG == SwGPrime)

        # Compare addition in different forms
        SwP = SWPoint.generator().scalar_mul(random_multiple)
        assert(SwP.is_on_curve())
        assert((SwP + SwG).to_montgomery() == SwG.to_montgomery() + SwP.to_montgomery())
        assert((SwP + SwG).to_montgomery().to_twisted_edwards() == SwG.to_montgomery().to_twisted_edwards() + SwP.to_montgomery().to_twisted_edwards())

        MontP = MontPoint.generator().scalar_mul(random_multiple)
        assert(MontP.is_on_curve())
        assert((MontP + MontG).to_short_weierstrass() == MontG.to_short_weierstrass() + MontP.to_short_weierstrass())
        assert((MontP + MontG).to_twisted_edwards() == MontG.to_twisted_edwards() + MontP.to_twisted_edwards())

        TwEdP = TwEdPoint.generator().scalar_mul(random_multiple)
        assert(TwEdP.is_on_curve())
        assert((TwEdP + TwEdG).to_montgomery() == TwEdG.to_montgomery() + TwEdP.to_montgomery())
        assert((TwEdP + TwEdG).to_montgomery().to_short_weierstrass() == TwEdG.to_montgomery().to_short_weierstrass() + TwEdP.to_montgomery().to_short_weierstrass())
        
        # Compare scaling in different forms
        TwEdGScaled = TwEdG.scalar_mul(scalar)
        assert(TwEdGScaled.is_on_curve())
        MontGScaled = MontG.scalar_mul(scalar)
        assert(MontGScaled.is_on_curve())
        SwGScaled = SwG.scalar_mul(scalar)
        assert(SwGScaled.is_on_curve())

        assert(SwGScaled.to_montgomery() == MontGScaled)
        assert(MontGScaled.to_short_weierstrass() == SwGScaled)
        assert(TwEdGScaled.to_montgomery() == MontGScaled)
        assert(MontGScaled.to_twisted_edwards() == TwEdGScaled)
        assert(TwEdGScaled.to_montgomery().to_short_weierstrass() == SwGScaled)
        assert(SwGScaled.to_montgomery().to_twisted_edwards() == TwEdGScaled)

if __name__ == '__main__':
    main()