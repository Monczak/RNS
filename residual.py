import math
import itertools
import functools


class Residual:
    def __init__(self, *args):
        if len(args) == 0:
            self.residual = None
            self.base = None
        else:
            if isinstance(args[0], tuple):
                self.residual = args[0]
            else:
                self.residual = self.to_residue(args[0], args[1])
            self.base = args[1]

    def __str__(self):
        return f"{' '.join(map(str, self.residual))} {str(self.base)}"

    @classmethod
    def is_valid_base(cls, base):
        return all([math.gcd(a, b) == 1 for a, b in itertools.combinations(base, 2)])

    @classmethod
    def to_residue(cls, number, base):
        if not cls.is_valid_base(base):
            raise Exception("Base is not valid")
        return tuple([number % b for b in base])

    def to_number(self, **kwargs):
        if "base" not in kwargs:
            base = self.base
        else:
            base = kwargs["base"]

        c = {}

        for i in range(2, len(base) + 1):
            c[i] = 1
            for j in range(1, i):
                u = self.inverse(base[j - 1], base[i - 1])
                c[i] = (c[i] * u) % base[i - 1]
        x = self.residual[0]
        for i in range(2, len(base) + 1):
            u = ((self.residual[i - 1] - x) * c[i]) % base[i - 1]
            x += u * functools.reduce(lambda a, b: a * b, base[:i - 1])
        return x

    @classmethod
    def inverse(cls, a, m):
        (b, x, y, n) = (m, 1, 0, 0)
        while a != 0:
            n = int(b / a)
            (a, b, x, y) = (b - n * a, a, y - n * x, x)
        return y % m

    @classmethod
    def sum_multiple(cls, residuals, base):
        result = Residual(0, base)
        for x in residuals:
            result += x
        return result

    @classmethod
    def to_mixed_radix(cls, residual, base, result=[]):
        if residual:
            a1 = residual[0]
            result.append(a1)

            a_i = [(a - a1) for a in residual[1:]]
            inverses = [cls.inverse(base[0], m) for m in base[1:]]
            a_j = [(a * m) % mi for (a, m, mi) in zip(a_i, inverses, base[1:])]
            return cls.to_mixed_radix(a_j, base[1:], result)
        else:
            result.reverse()
            return result

    @staticmethod
    def comparer(a, b):
        if a == b:
            return "="
        a_radix = Residual.to_mixed_radix(a.residual, a.base, [])
        b_radix = Residual.to_mixed_radix(b.residual, b.base, [])

        for i in range(len(a_radix)):
            if a_radix[i] > b_radix[i]:
                return ">"
            elif a_radix[i] < b_radix[i]:
                return "<"

    def __add__(self, other):
        result = Residual()

        if not isinstance(other, Residual):
            other = Residual(other, self.base)

        if self.base != other.base:
            raise Exception("Numbers must have same base")

        result.residual = tuple([(self.residual[i] + other.residual[i]) % self.base[i] for i in range(len(self.base))])
        result.base = self.base
        return result

    def __sub__(self, other):
        result = Residual()

        if not isinstance(other, Residual):
            other = Residual(other, self.base)

        if self.base != other.base:
            raise Exception("Numbers must have same base")

        result.residual = tuple([(self.residual[i] - other.residual[i]) % self.base[i] for i in range(len(self.base))])
        result.base = self.base
        return result

    def __mul__(self, other):
        result = Residual()

        if not isinstance(other, Residual):
            other = Residual(other, self.base)

        if self.base != other.base:
            raise Exception("Numbers must have same base")

        result.residual = tuple([(self.residual[i] * other.residual[i]) % self.base[i] for i in range(len(self.base))])
        result.base = self.base
        return result

    def __floordiv__(self, other):
        if self == other:
            return Residual(1, self.base)
        x = []
        u = []
        n = Residual(1, self.base)
        while self >= n * other:
            x.append(n)
            n *= 2
        x.reverse()
        for xn in x:
            if (self.sum_multiple(u, self.base) + xn) * other <= self:
                u.append(xn)
        return self.sum_multiple(u, self.base)

    def __lt__(self, other):
        if isinstance(other, Residual):
            return self.comparer(self, other) == "<"
        else:
            return self.comparer(self, Residual(other, self.base)) == "<"

    def __gt__(self, other):
        if isinstance(other, Residual):
            return self.comparer(self, other) == ">"
        else:
            return self.comparer(self, Residual(other, self.base)) == ">"

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other

    def __eq__(self, other):
        if isinstance(other, Residual):
            return self.residual == other.residual and self.base == other.base
        else:
            return self.residual == Residual(other, self.base).residual
