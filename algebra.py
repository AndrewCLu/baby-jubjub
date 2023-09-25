# Some basic finite field algebra
class FElt:
    def __init__(self, field, value):
        if not isinstance(field, GF):
            raise ValueError('Field must be an instance of GaloisField')
        if value >= field.order or value < 0:
            raise ValueError('Value must be between 0 and the field order')
        self.field = field
        self.value = value
    
    def __add__(self, other):
        if self.field.order != other.field.order:
            raise TypeError('Cannot add elements from different fields')
        return self.field((self.value + other.value) % self.field.order)
    
    def __sub__(self, other):
        if self.field.order != other.field.order:
            raise TypeError('Cannot subtract elements from different fields')
        return self.field((self.value - other.value) % self.field.order)
    
    def __mul__(self, other):
        if self.field.order != other.field.order:
            raise TypeError('Cannot multiply elements from different fields')
        return self.field((self.value * other.value) % self.field.order)
    
    def __truediv__(self, other):
        if self.field.order != other.field.order:
            raise TypeError('Cannot divide elements from different fields')
        return self * self.field(pow(other.value, self.field.order - 2, self.field.order))

    def __neg__(self):
        if self.value == 0:
            return self
        return self.field(self.field.order - self.value)
    
    def __eq__(self, other):
        if self.field.order != other.field.order:
            raise TypeError('Cannot compare elements from different fields')
        return self.value == other.value
    
    def __str__(self):
        return str(self.value)

class GF:
    def __init__(self, order):
        self.order = order
    
    def __call__(self, value):
        return FElt(self, value)