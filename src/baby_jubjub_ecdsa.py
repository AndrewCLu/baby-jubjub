# Implementation of ECDSA over Baby Jubjub
from algebra import GF
from baby_jubjub import BabyJubjubPoint, SWPoint, MontPoint, TwEdPoint
from ecdsa.numbertheory import SquareRootError

# Note we use the prime order subgroup generated by the base point of Baby Jubjub
cofactor = 8
order = BabyJubjubPoint.order // cofactor
Fr = BabyJubjubPoint.Fr
Fn = GF(order)

def keygen(representation, seed):
    assert(representation in [SWPoint, MontPoint, TwEdPoint])

    # Generate a private key
    if seed <= 0 or seed >= order:
        raise ValueError("Seed must be in the range [1, order - 1]")
    priv = Fn(seed).value
    
    # Generate the public key
    pub = representation.base().scalar_mul(priv)
    
    return (priv, pub)

def sign(digest, privKey, representation):
    assert(representation in [SWPoint, MontPoint, TwEdPoint])
    assert(isinstance(digest, int))
    assert(isinstance(privKey, int))
    assert(privKey >= 0 and privKey < order)

    # Generate a random nonce
    k = 10
    
    # Generate the random point
    R = representation.base().scalar_mul(k)
    r = Fn(R.x.value)
    if r == Fn(0):
        raise ValueError("Failed to generate a valid signature. Try again with a different nonce.")
    
    # Compute the signature
    s = (Fn(digest) + r * Fn(privKey)) / Fn(k)
    if s == Fn(0):
        raise ValueError("Failed to generate a valid signature. Try again with a different nonce.")
    
    return (r.value, s.value)

# Standard ECDSA verification algorithm
def verify(digest, pubKey, r, s, representation):
    assert(representation in [SWPoint, MontPoint, TwEdPoint])
    assert(isinstance(digest, int))
    assert(isinstance(pubKey, representation))
    assert(isinstance(r, int))
    assert(isinstance(s, int))

    if r <= 0 or r >= order or s <= 0 or s >= order:
        return False

    u_1 = Fn(digest) / Fn(s)
    u_2 = Fn(r) / Fn(s)

    pt = representation.base().scalar_mul(u_1.value) + pubKey.scalar_mul(u_2.value)
    if pt.is_infinity():
        return False
    
    return Fn(r) == Fn(pt.x.value)

# Verify a signature taking in R as advice
# R is the point computed as R = k * G
# This allows us to verify a signature in a different representation than the one it was signed in
# R can be recovered from r in the original coordinate system, and then transformed to the new coordinate system
# The verification equation is taken from Efficient ECDSA: https://personaelabs.org/posts/efficient-ecdsa-1/
def verify_with_advice(digest, pubKey, r, s, R, representation):
    assert(representation in [SWPoint, MontPoint, TwEdPoint])
    assert(isinstance(digest, int))
    assert(isinstance(pubKey, representation))
    assert(isinstance(r, int))
    assert(isinstance(s, int))
    assert(isinstance(R, representation))

    if r <= 0 or r >= order or s <= 0 or s >= order:
        return False

    sR = R.scalar_mul(s)
    mG = representation.base().scalar_mul(digest)
    rQa = pubKey.scalar_mul(r)
    return sR == mG + rQa

def recover_public_key(digest, r, s, representation):
    if r <= 0 or r >= order or s <= 0 or s >= order:
        return False
    
    possible_points = representation.recover_from_x(r)
    possible_pub_keys = []
    for pt in possible_points:
        u_1 = Fn(0) - Fn(digest) / Fn(r)
        u_2 = Fn(s) / Fn(r)

        pub_key = representation.base().scalar_mul(u_1.value) + pt.scalar_mul(u_2.value)
        if (verify(digest, pub_key, r, s, representation)):
            possible_pub_keys.append(pub_key)

    return possible_pub_keys