from baby_jubjub_ecdsa import keygen, sign, verify
from baby_jubjub import BabyJubjubPoint, SWPoint, TwEdPoint
from ecdsa import ellipticcurve, curves
from ecdsa.ellipticcurve import PointEdwards
from ecdsa.ecdsa import Public_key, Private_key

# Verifies our implementation of ECDSA over Baby Jubjub matches python-ecdsa
# https://github.com/tlsfuzzer/python-ecdsa/tree/master
def main():
    # Short Weierstrass form
    p = BabyJubjubPoint.p
    prime_subgroup_order = BabyJubjubPoint.prime_subgroup_order

    SWellipticCurve = ellipticcurve.CurveFp(p, SWPoint.a.value, SWPoint.b.value)
    TwEdellipticCurve = ellipticcurve.CurveEdTw(p, TwEdPoint.A.value, TwEdPoint.d.value)

    def getTwEdPointFromxy(x, y):
        return PointEdwards(
            TwEdellipticCurve, x, y, 1, x * y, prime_subgroup_order, False
        )
    
    SWB = ellipticcurve.Point(SWellipticCurve, SWPoint.base().x.value, SWPoint.base().y.value, order=prime_subgroup_order)
    SWcurve = curves.Curve('babyjubjub-sw', SWellipticCurve, SWB, None)
    TwEdB = getTwEdPointFromxy(TwEdPoint.base().x.value, TwEdPoint.base().y.value)
    TwEdcurve = curves.Curve('babyjubjub-twed', TwEdellipticCurve, TwEdB, None)

    # Test 1: Verify conversion between representations
    mySWB = SWPoint.base()
    assert(ellipticcurve.Point(SWellipticCurve, mySWB.x.value, mySWB.y.value, order=prime_subgroup_order) == SWB)
    myTwEdB = TwEdPoint.base()
    assert(getTwEdPointFromxy(myTwEdB.x.value, myTwEdB.y.value) == TwEdB)

    mySWB = SWPoint(SWB.x(), SWB.y())
    myTwEdB = mySWB.to_montgomery().to_twisted_edwards()
    assert(getTwEdPointFromxy(myTwEdB.x.value, myTwEdB.y.value) == TwEdB)

    myTwEdB = TwEdPoint(TwEdB.x(), TwEdB.y())
    mySWB = myTwEdB.to_montgomery().to_short_weierstrass()
    assert(ellipticcurve.Point(SWellipticCurve, mySWB.x.value, mySWB.y.value, order=prime_subgroup_order) == SWB)

    # Test 2: Verify operations and conversions
    scalar = 1023
    myScaledSWB = SWPoint(SWB.x(), SWB.y()).scalar_mul(scalar)
    myScaledTwEdB = myScaledSWB.to_montgomery().to_twisted_edwards()
    assert(getTwEdPointFromxy(myScaledTwEdB.x.value, myScaledTwEdB.y.value) == TwEdB * scalar)

    scalar = 1234
    myScaledTwEdB = TwEdPoint(TwEdB.x(), TwEdB.y()).scalar_mul(scalar)
    myScaledSWB = myScaledTwEdB.to_montgomery().to_short_weierstrass()
    assert(ellipticcurve.Point(SWellipticCurve, myScaledSWB.x.value, myScaledSWB.y.value, order=prime_subgroup_order) == SWB * scalar)

    # Test 3: Verify signatures are the same in SW form
    secrets = [100, 200, 300]
    digests = [1000, 2000, 3000]
    random_ks = [10000, 20000, 30000]

    for secret in secrets:
        for digest in digests:
            for random_k in random_ks:
                # print(f"Verifying same SW signature is generated for secret {secret}, digest {digest}, random_k {random_k}")

                pubKey = Public_key(SWB, SWB * secret)
                privKey = Private_key(pubKey, secret)
                (myPrivKey, myPubKey) = keygen(SWPoint, secret)
                assert(myPrivKey == privKey.secret_multiplier)
                assert(myPubKey.x.value == pubKey.point.x())
                assert(myPubKey.y.value == pubKey.point.y())

                sig = privKey.sign(digest, random_k)
                r = sig.r
                s = sig.s
                (myR, myS) = sign(SWPoint, digest, myPrivKey, random_k)
                assert(r == myR)
                assert(s == myS)

                assert(pubKey.verifies(digest, sig))
                assert(verify(SWPoint, digest, myPubKey, r, s))

if __name__ == '__main__':
    main()