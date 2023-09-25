from algebra import GF

class BabyJubjubPoint:
    # Base field
    p=21888242871839275222246405745257275088548364400416034343698204186575808495617
    Fr=GF(p)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError("Can only compare points in the same representation.")
        if self.is_infinity() and other.is_infinity():
            return True
        if self.is_infinity() or other.is_infinity():
            return False
        return self.x == other.x and self.y == other.y

    def scalar_mul(self, scalar):
        if not isinstance(scalar, int):
            raise TypeError("Scalar must be an integer.")
        if scalar < 0:
            raise ValueError("Scalar must be non-negative.")
        if scalar == 0:
            return self.__class__.point_at_infinity()
        if scalar == 1:
            return self
        if scalar % 2 == 0:
            return (self + self).scalar_mul(scalar // 2)
        else:
            return self + self.scalar_mul(scalar - 1)
        
# A Baby Jubjub point represented in Short Weierstrass form
class SWPoint(BabyJubjubPoint):
    # Short Weierstrass parameters
    a=BabyJubjubPoint.Fr(7296080957279758407415468581752425029516121466805344781232734728849116493472)
    b=BabyJubjubPoint.Fr(16213513238399463127589930181672055621146936592900766180517188641980520820846)

    def __init__(self, x, y):
        # Point at infinity
        if x is None and y is None:  
            self.x = self.y = None
            return

        if isinstance(x, int):
            self.x = self.Fr(x)
        else:
            self.x = x
        if isinstance(y, int):
            self.y = self.Fr(y)
        else:
            self.y = y

        if not self.is_on_curve():
            raise ValueError(f"The point ({x}, {y}) is not on the curve.")
    
    def __add__(self, other):
        if not isinstance(other, SWPoint):
            raise TypeError("Can only add Short Weierstrass points to other Short Weierstrass points.")
        
        if self.is_infinity():
            return other
        
        if other.is_infinity():
            return self
        
        # Adding a point and its inverse gives infinity
        if self.x == other.x and self.y == -other.y:
            return SWPoint(None, None)
        
        x1 = self.x
        y1 = self.y
        x2 = other.x
        y2 = other.y
        
        if x1 == x2 and y1 == y2:
            lam = (self.Fr(3) * x1 * x1 + self.a) / (self.Fr(2) * y1)
        else:
            lam = (y2 - y1) / (x2 - x1)
        
        x3 = lam * lam - x1 - x2
        y3 = lam * (x1 - x3) - y1

        return SWPoint(x3, y3)
    
    def __str__(self):
        if self.is_infinity():
            return "Inf"
        return f"SW: ({self.x}, {self.y})"
    
    def is_infinity(self):
        return self.x is None and self.y is None
    
    def point_at_infinity():
        return SWPoint(None, None)

    def is_on_curve(self):
        if self.is_infinity():
            return True
        
        lhs = self.y * self.y
        rhs = self.x * self.x * self.x + self.a * self.x + self.b

        return lhs == rhs
    
    def to_montgomery(self):
        if self.is_infinity():
            return MontPoint.point_at_infinity()
        
        alpha = MontPoint.alpha
        beta = MontPoint.beta
        nx = (self.x - alpha) / beta
        ny = self.y / beta
        return MontPoint(nx, ny)
    
# A Baby Jubjub point represented in Montgomery form
class MontPoint(BabyJubjubPoint):
    # Montgomery parameters
    MontA = BabyJubjubPoint.Fr(168698)
    MontB = BabyJubjubPoint.Fr(1)
    alpha = MontA / BabyJubjubPoint.Fr(3)
    beta = BabyJubjubPoint.Fr(1) / MontB

    def __init__(self, x, y):
        # Point at infinity
        if x is None and y is None:  
            self.x = self.y = None
            return

        if isinstance(x, int):
            self.x = self.Fr(x)
        else:
            self.x = x
        if isinstance(y, int):
            self.y = self.Fr(y)
        else:
            self.y = y

        if not self.is_on_curve():
            raise ValueError(f"The point ({x}, {y}) is not on the curve.")
    
    def __add__(self, other):
        if not isinstance(other, MontPoint):
            raise TypeError("Can only add Montgomery points to other Montgomery points.")
        
        if self.is_infinity():
            return other
        
        if other.is_infinity():
            return self
        
        # Adding a point and its inverse gives infinity
        if self.x == other.x and self.y == -other.y:
            return MontPoint(None, None)
        
        x1 = self.x
        y1 = self.y
        x2 = other.x
        y2 = other.y
        
        if x1 == x2 and y1 == y2:
            lam = (self.Fr(3) * x1 * x1 + self.Fr(2) * self.MontA * x1 + self.Fr(1)) / (self.Fr(2) * self.MontB * y1)
        else:
            lam = (y2 - y1) / (x2 - x1)
        
        x3 = self.MontB * lam * lam - self.MontA - x1 - x2
        y3 = (self.Fr(2) * x1 + x2 + self.MontA) * lam - self.MontB * lam * lam * lam - y1

        return MontPoint(x3, y3)
    
    def __str__(self):
        if self.is_infinity():
            return "Inf"
        return f"Mont: ({self.x}, {self.y})"
    
    def point_at_infinity():
        return MontPoint(None, None)
    
    def is_infinity(self):
        return self.x is None and self.y is None
    
    def is_on_curve(self):
        if self.is_infinity():
            return True
        
        lhs = self.MontB * self.y * self.y
        rhs = self.x * self.x * self.x + self.MontA * self.x * self.x + self.x

        return lhs == rhs
    
    def to_short_weierstrass(self):
        if self.is_infinity():
            return SWPoint.point_at_infinity()
        
        nx = (self.x + self.MontA / self.Fr(3)) / self.MontB
        ny = self.x / self.MontB
        return SWPoint(nx, ny)
    
    def to_twisted_edwards(self):
        if self.is_infinity():
            return TwEdPoint.point_at_infinity()
        
        nx = self.x / self.y
        ny = (self.x - self.Fr(1)) / (self.x + self.Fr(1))
        return TwEdPoint(nx, ny)
    
# A Baby Jubjub point represented in Twisted Edwards form
class TwEdPoint(BabyJubjubPoint):
    # Twisted Edwards parameters
    TwEdA = BabyJubjubPoint.Fr(168700)
    TwEdd = BabyJubjubPoint.Fr(168696)

    def __init__(self, x, y):
        # No special point at infinity in Twisted Edwards form
        if isinstance(x, int):
            self.x = self.Fr(x)
        else:
            self.x = x
        if isinstance(y, int):
            self.y = self.Fr(y)
        else:
            self.y = y

        if not self.is_on_curve():
            raise ValueError(f"The point ({x}, {y}) is not on the curve.")
    
    def __add__(self, other):
        if not isinstance(other, TwEdPoint):
            raise TypeError("Can only add Twisted Edwards points to other Twisted Edwards points.")
        
        x1 = self.x
        y1 = self.y
        x2 = other.x
        y2 = other.y
        
        x3 = (x1 * y2 + y1 * x2) / (self.Fr(1) + self.TwEdd * x1 * x2 * y1 * y2)
        y3 = (y1 * y2 - self.TwEdA * x1 * x2) / (self.Fr(1) - self.TwEdd * x1 * x2 * y1 * y2)

        return TwEdPoint(x3, y3)
    
    def __str__(self):
        return f"TwEd: ({self.x}, {self.y})"
    
    def point_at_infinity():
        return TwEdPoint(BabyJubjubPoint.Fr(0), BabyJubjubPoint.Fr(1))
    
    def is_infinity(self):
        return self.x == self.Fr(0) and self.y == self.Fr(1)

    def is_on_curve(self):
        lhs = self.TwEdA * self.x * self.x + self.y * self.y
        rhs = self.Fr(1) + self.TwEdd * self.x * self.x * self.y * self.y

        return lhs == rhs
    
    def to_montgomery(self):
        if self.is_infinity():
            return MontPoint.point_at_infinity()
        
        nx = (self.Fr(1) + self.y) / (self.Fr(1) - self.y)
        ny = (self.Fr(1) + self.y) / ((self.Fr(1) - self.y) * self.x)
        return MontPoint(nx, ny)