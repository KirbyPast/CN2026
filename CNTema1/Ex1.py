m = 1
u = pow(10,-m)

while 1.0 + u != 1.0:
    m = m+1
    u = pow(10,-m)

print(m)
