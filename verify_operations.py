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
    alpha = MontA / Fr(3)
    beta = Fr(1) / MontB

    # Test some conversions
    SwGx = 14414009007687342025526645003307639786191886886413750648631138442071909631647
    SwGy = 14577268218881899420966779687690205425227431577728659819975198491127179315626
    for i in range(1000):
        SwG = SWPoint(SwGx, SwGy).scalar_mul(i)
        print(i)
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
        
        # Operations in Twisted Edwards form
        TwEdGScaled = TwEdG.scalar_mul(11)
        assert(TwEdGScaled.is_on_curve())
        MontGScaled = MontG.scalar_mul(11)
        assert(MontGScaled.is_on_curve())
        assert(TwEdGScaled.to_montgomery() == MontGScaled)

        # Convert from Montgomery to Short Weierstrass
        # try:
        #     SwGPrime = MontGPrime.to_short_weierstrass()
        #     assert(SwGPrime.is_on_curve())
        #     assert(SwG == SwGPrime)
        #     # print("NO ERROR")
        # except:
        #     # print("error")
        #     pass
    

if __name__ == '__main__':
    main()