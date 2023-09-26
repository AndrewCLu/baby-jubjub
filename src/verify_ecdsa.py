from baby_jubjub_ecdsa import keygen, sign, verify, recover_public_key, verify_with_advice
from baby_jubjub import SWPoint, MontPoint, TwEdPoint
import json

def main():
    print("Verifying ECDSA over Baby Jubjub...")
    # Generate some signatures
    # signatures = []
    # for seed in range(100, 301, 100):
    #     for digest in range(1000, 1301, 100):
    #         for representation in [SWPoint, MontPoint, TwEdPoint]:
    #             print(f"Testing {representation.__name__} with seed {seed} and digest {digest}")
    #             (privKey, pubKey) = keygen(representation, seed)
    #             (r, s) = sign(digest, privKey, representation)
    #             assert(verify(digest, pubKey, r, s, representation))
    #             sig = [seed, digest, representation.__name__, privKey, pubKey.x.value, pubKey.y.value, r, s]
    #             signatures.append(sig)

    # with open('signatures.json', 'w') as f:
    #     json.dump(signatures, f)

    with open('signatures.json', 'r') as f:
        signatures = json.load(f)

    for sig in signatures:
        seed, digest, representationName, privKey, pubKeyX, pubKeyY, r, s = sig
        representationMap = {
            'SWPoint': SWPoint,
            'MontPoint': MontPoint,
            'TwEdPoint': TwEdPoint
        }
        representation = representationMap[representationName]
        pubKey = representation(pubKeyX, pubKeyY)
        print(f"Verifying {representation.__name__} with seed {seed} and digest {digest}")
        if representation.__name__ == 'SWPoint':
            print("Verify SWPoint in Short Weierstrass form: ", verify(digest, pubKey, r, s, SWPoint))

            verifiedMont = False
            verifiedTwEd = False
            pubKeyMont = pubKey.to_montgomery()
            pubKeyTwEd = pubKeyMont.to_twisted_edwards()
            Rs = representation.recover_from_x(r)
            for R in Rs:
                RMont = R.to_montgomery()
                RTwed = RMont.to_twisted_edwards()
                if verify_with_advice(digest, pubKeyMont, r, s, RMont, MontPoint):
                    verifiedMont = True
                if verify_with_advice(digest, pubKeyTwEd, r, s, RTwed, TwEdPoint):
                    verifiedTwEd = True

            print("Verify SWPoint in Montgomery form: ", verifiedMont)
            print("Verify SWPoint in Twisted Edwards form: ", verifiedTwEd)
        elif representation.__name__ == 'MontPoint':
            print("Verify MontPoint in Montgomery form: ", verify(digest, pubKey, r, s, MontPoint))

            verifiedSW = False
            verifiedTwEd = False
            pubKeySW = pubKey.to_short_weierstrass()
            pubKeyTwEd = pubKey.to_twisted_edwards()
            Rs = representation.recover_from_x(r)
            for R in Rs:
                RSW = R.to_short_weierstrass()
                RTwed = R.to_twisted_edwards()
                if verify_with_advice(digest, pubKeySW, r, s, RSW, SWPoint):
                    verifiedSW = True
                if verify_with_advice(digest, pubKeyTwEd, r, s, RTwed, TwEdPoint):
                    verifiedTwEd = True

            print("Verify MontPoint in Short Weierstrass form: ", verifiedSW)
            print("Verify MontPoint in Twisted Edwards form: ", verifiedTwEd)
        elif representation.__name__ == 'TwEdPoint':
            print("Verify TwEdPoint in Twisted Edwards form: ", verify(digest, pubKey, r, s, TwEdPoint))
            
            verifiedSW = False
            verifiedMont = False
            pubKeyMont = pubKey.to_montgomery()
            pubKeySW = pubKeyMont.to_short_weierstrass()
            Rs = representation.recover_from_x(r)
            for R in Rs:
                RMont = R.to_montgomery()
                RSW = RMont.to_short_weierstrass()
                if verify_with_advice(digest, pubKeySW, r, s, RSW, SWPoint):
                    verifiedSW = True
                if verify_with_advice(digest, pubKeyMont, r, s, RMont, MontPoint):
                    verifiedMont = True

            print("Verify TwEdPoint in Short Weierstrass form: ", verifiedSW)
            print("Verify TwEdPoint in Montgomery form: ", verifiedMont)
        # assert(verify(digest, pubKey, r, s, representation))
        # assert(pubKey in recover_public_key(digest, r, s, representation))
        print(f"Verified {representation.__name__} with seed {seed} and digest {digest}\n")

    print("ECDSA over Baby Jubjub verified!")

if __name__ == "__main__":
    main()