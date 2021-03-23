from residual import Residual

base = (31, 32, 33)

a = Residual(23, base)
b = Residual(18, base)
c = Residual(23, base)

print("a + b", a + b)
print("a - b", a - b)
print("a * b", a * b)
print("a // b", a // b)

print(a.to_number(), b.to_number(), (a + b).to_number())

print("a > b is", a > b, "a < b is", a < b)
print(a > c, a < c, a >= c, a <= c, a == c)
print(Residual.to_mixed_radix(Residual(320, base).residual, base))
print((Residual(3200, base) // Residual(8, base)).to_number())
