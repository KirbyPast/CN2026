from random import random

x = 1.0
u = pow(10,-15)
y = u / 10
z = u / 10
if (x+y) + z == x + (y+z):
    print("Adunarea calculator este asociativa")
else:
    print("Adunarea calculator NU este asociativa")
counter = 0
while (x * y) * z == x * (y * z):
    counter = counter + 1
    x=random()
    y=random()
    z=random()

print(x.__str__() + ", " + y.__str__() + ", " + z.__str__() + " (" + counter.__str__() + " incercari)")