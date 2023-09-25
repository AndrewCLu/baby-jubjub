from baby_jubjub import SWPoint, MontPoint, TwEdPoint

# Verifies different forms of the Baby Jubjub curve are equivalent
# Reference: https://www-fourier.univ-grenoble-alpes.fr/mphell/doc-v5/conversion_weierstrass_edwards.html
def main():
    print("Verifying operations on the Baby Jubjub curve...")
    # Test some operations using Short Weierstrass Generator point
    scalar = 1023
    random_multiple = 1234567890
    num_points_to_test = 100
    for i in range(0, SWPoint.order, SWPoint.order // num_points_to_test):
        print(f"Verifying conversions for {i}G")
        # Generate a point on the curve
        SWG = SWPoint.generator().scalar_mul(i)
        assert(SWG.is_on_curve())

        # Convert from Short Weierstrass to Montgomery
        MontG = SWG.to_montgomery()
        assert(MontG.is_on_curve())

        # Convert from Montgomery to Twisted Edwards
        TwEdG = MontG.to_twisted_edwards()
        assert(TwEdG.is_on_curve())

        # Convert from Twisted Edwards to Montgomery
        MontGPrime = TwEdG.to_montgomery()
        assert(MontGPrime.is_on_curve())
        assert(MontG == MontGPrime)

        # Convert from Montgomery to Short Weierstrass
        SWGPrime = MontGPrime.to_short_weierstrass()
        assert(SWGPrime.is_on_curve())
        assert(SWG == SWGPrime)

        # Compare addition in different forms
        SWP = SWPoint.generator().scalar_mul(random_multiple)
        assert(SWP.is_on_curve())
        assert((SWP + SWG).to_montgomery() == SWG.to_montgomery() + SWP.to_montgomery())
        assert((SWP + SWG).to_montgomery().to_twisted_edwards() == SWG.to_montgomery().to_twisted_edwards() + SWP.to_montgomery().to_twisted_edwards())

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
        SWGScaled = SWG.scalar_mul(scalar)
        assert(SWGScaled.is_on_curve())

        assert(SWGScaled.to_montgomery() == MontGScaled)
        assert(MontGScaled.to_short_weierstrass() == SWGScaled)
        assert(TwEdGScaled.to_montgomery() == MontGScaled)
        assert(MontGScaled.to_twisted_edwards() == TwEdGScaled)
        assert(TwEdGScaled.to_montgomery().to_short_weierstrass() == SWGScaled)
        assert(SWGScaled.to_montgomery().to_twisted_edwards() == TwEdGScaled)
    
    print("Curve operations verified!")

if __name__ == '__main__':
    main()