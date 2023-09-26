from algebra import GF
from ecdsa.numbertheory import SquareRootError

class BabyJubjubPoint:
    # Base field
    p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
    Fr = GF(p)
    order = 21888242871839275222246405745257275088614511777268538073601725287587578984328
    cofactor = 8
    prime_subgroup_order = order // cofactor

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
            return self.__class__.infinity()
        if scalar == 1:
            return self
        if scalar % 2 == 0:
            return (self + self).scalar_mul(scalar // 2)
        else:
            return self + self.scalar_mul(scalar - 1)
        
# A Baby Jubjub point represented in Short Weierstrass form
class SWPoint(BabyJubjubPoint):
    Fr = BabyJubjubPoint.Fr
    # Short Weierstrass parameters
    a = Fr(7296080957279758407415468581752425029516121466805344781232734728849116493472)
    b = Fr(16213513238399463127589930181672055621146936592900766180517188641980520820846)
    Gx = Fr(7296080957279758407415468581752425029516121466805344781232734728858602888112)
    Gy = Fr(4258727773875940690362607550498304598101071202821725296872974770776423442226)
    Bx = Fr(14414009007687342025526645003307639786191886886413750648631138442071909631647)
    By = Fr(14577268218881899420966779687690205425227431577728659819975198491127179315626)

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
    
    def infinity():
        return SWPoint(None, None)
    
    def generator():
        return SWPoint(SWPoint.Gx, SWPoint.Gy)
    
    def base():
        return SWPoint(SWPoint.Bx, SWPoint.By)
    
    def recover_from_x(x_int):
        if not isinstance(x_int, int):
            raise TypeError("x must be an integer.")
        if x_int < 0 or x_int >= SWPoint.p:
            raise ValueError("x must be in the field.")
        
        possible_points = []
        # r is modded by the order of the subgroup, so we must try all possible values of r
        for m in range(SWPoint.cofactor):
            x = SWPoint.Fr(x_int) + SWPoint.Fr(m) * SWPoint.Fr(SWPoint.prime_subgroup_order)
            y2 = x * x * x + SWPoint.a * x + SWPoint.b
            try:
                y = y2.sqrt()
                possible_points += [SWPoint(x, y), SWPoint(x, -y)]
            except SquareRootError:
                pass
        
        return possible_points

    def is_on_curve(self):
        if self.is_infinity():
            return True
        
        lhs = self.y * self.y
        rhs = self.x * self.x * self.x + self.a * self.x + self.b

        return lhs == rhs
    
    def to_montgomery(self):
        if self.is_infinity():
            return MontPoint.infinity()
        
        nx = self.x - MontPoint.alpha
        ny = self.y
        return MontPoint(nx, ny)
    
# A Baby Jubjub point represented in Montgomery form
class MontPoint(BabyJubjubPoint):
    Fr = BabyJubjubPoint.Fr
    # Montgomery parameters
    A = Fr(168698)
    B = Fr(1)
    alpha = A / Fr(3)
    Gx = Fr(7)
    Gy = Fr(4258727773875940690362607550498304598101071202821725296872974770776423442226)
    Bx = Fr(7117928050407583618111176421555214756675765419608405867398403713213306743542)
    By = Fr(14577268218881899420966779687690205425227431577728659819975198491127179315626)

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
            lam = (self.Fr(3) * x1 * x1 + self.Fr(2) * self.A * x1 + self.Fr(1)) / (self.Fr(2) * self.B * y1)
        else:
            lam = (y2 - y1) / (x2 - x1)
        
        x3 = self.B * lam * lam - self.A - x1 - x2
        y3 = (self.Fr(2) * x1 + x2 + self.A) * lam - self.B * lam * lam * lam - y1

        return MontPoint(x3, y3)
    
    def __str__(self):
        if self.is_infinity():
            return "Inf"
        return f"Mont: ({self.x}, {self.y})"
    
    def is_infinity(self):
        return self.x is None and self.y is None
    
    def infinity():
        return MontPoint(None, None)
    
    def generator():
        return MontPoint(MontPoint.Gx, MontPoint.Gy)
    
    def base():
        return MontPoint(MontPoint.Bx, MontPoint.By)
    
    def is_on_curve(self):
        if self.is_infinity():
            return True
        
        lhs = self.B * self.y * self.y
        rhs = self.x * self.x * self.x + self.A * self.x * self.x + self.x

        return lhs == rhs
        
    def recover_from_x(x_int):
        if not isinstance(x_int, int):
            raise TypeError("x must be an integer.")
        if x_int < 0 or x_int >= MontPoint.p:
            raise ValueError("x must be in the field.")
        
        possible_points = []
        # r is modded by the order of the subgroup, so we must try all possible values of r
        for m in range(MontPoint.cofactor):
            x = MontPoint.Fr(x_int) + MontPoint.Fr(m) * MontPoint.Fr(MontPoint.prime_subgroup_order)
            y2 = (x * x * x + MontPoint.A * x * x + x) / MontPoint.B
            try:
                y = y2.sqrt()
                possible_points += [MontPoint(x, y), MontPoint(x, -y)]
            except SquareRootError:
                pass
            
        return possible_points
    
    def to_short_weierstrass(self):
        if self.is_infinity():
            return SWPoint.infinity()
        
        nx = (self.x + self.A / self.Fr(3)) / self.B
        ny = self.y / self.B
        return SWPoint(nx, ny)
    
    def to_twisted_edwards(self):
        if self.is_infinity():
            return TwEdPoint.infinity()
        
        nx = self.x / self.y
        ny = (self.x - self.Fr(1)) / (self.x + self.Fr(1))
        return TwEdPoint(nx, ny)
    
# A Baby Jubjub point represented in Twisted Edwards form
class TwEdPoint(BabyJubjubPoint):
    Fr = BabyJubjubPoint.Fr
    # Twisted Edwards parameters
    A = Fr(168700)
    d = Fr(168696)
    Gx = Fr(995203441582195749578291179787384436505546430278305826713579947235728471134)
    Gy = Fr(5472060717959818805561601436314318772137091100104008585924551046643952123905)
    Bx = Fr(5299619240641551281634865583518297030282874472190772894086521144482721001553)
    By = Fr(16950150798460657717958625567821834550301663161624707787222815936182638968203)

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
        
        x3 = (x1 * y2 + y1 * x2) / (self.Fr(1) + self.d * x1 * x2 * y1 * y2)
        y3 = (y1 * y2 - self.A * x1 * x2) / (self.Fr(1) - self.d * x1 * x2 * y1 * y2)

        return TwEdPoint(x3, y3)
    
    def __str__(self):
        return f"TwEd: ({self.x}, {self.y})"
    
    def is_infinity(self):
        return self.x == self.Fr(0) and self.y == self.Fr(1)
    
    def infinity():
        return TwEdPoint(BabyJubjubPoint.Fr(0), BabyJubjubPoint.Fr(1))
    
    def generator():
        return TwEdPoint(TwEdPoint.Gx, TwEdPoint.Gy)
    
    def base():
        return TwEdPoint(TwEdPoint.Bx, TwEdPoint.By)

    def is_on_curve(self):
        lhs = self.A * self.x * self.x + self.y * self.y
        rhs = self.Fr(1) + self.d * self.x * self.x * self.y * self.y

        return lhs == rhs
    
    def recover_from_x(x_int):
        if not isinstance(x_int, int):
            raise TypeError("x must be an integer.")
        if x_int < 0 or x_int >= TwEdPoint.p:
            raise ValueError("x must be in the field.")
        
        possible_points = []
        # r is modded by the order of the subgroup, so we must try all possible values of r
        for m in range(MontPoint.cofactor):
            x = TwEdPoint.Fr(x_int) + TwEdPoint.Fr(m) * TwEdPoint.Fr(TwEdPoint.prime_subgroup_order)
            y2 = (TwEdPoint.A * x * x - TwEdPoint.Fr(1)) / (TwEdPoint.d * x * x - TwEdPoint.Fr(1))
            try:
                y = y2.sqrt()
                possible_points += [TwEdPoint(x, y), TwEdPoint(x, -y)]
            except SquareRootError:
                pass
            
        return possible_points
    
    def to_montgomery(self):
        if self.is_infinity():
            return MontPoint.infinity()
        
        nx = (self.Fr(1) + self.y) / (self.Fr(1) - self.y)
        ny = (self.Fr(1) + self.y) / ((self.Fr(1) - self.y) * self.x)
        return MontPoint(nx, ny)