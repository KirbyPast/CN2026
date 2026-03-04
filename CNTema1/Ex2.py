from Ex1 import precizia_masina
from random import random

if __name__ == "__main__":
    x = 1.0
    u = precizia_masina()
    y = u / 10
    z = u / 10

    if (x + y) + z == x + (y + z):
        print("Adunarea calculator este asociativa")
    else:
        print("Adunarea calculator NU este asociativa")

    counter = 0

    while (x * y) * z == x * (y * z):
        counter = counter + 1
        x = random()
        y = random()
        z = random()

    print(f"{x}, {y}, {z} ({counter} incercari)")
