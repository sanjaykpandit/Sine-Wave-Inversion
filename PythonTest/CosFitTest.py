
import numpy as np
from CosFit import CosFit
import time

start = 0
stop = 1
n = 1000

a = 10.0
f = 1000; omega = 2.0 * np.pi * f;
p = np.pi / 4
d = 0.5 
t = []
x = []

m = 201
for i in range(m):
    t.append(np.linspace(start, stop, n))
    x.append(a * np.cos(omega*t[i]+p+0.01*i) + d)

st = time.time()

for i in range(m):
    cs = CosFit(t[i], x[i])
    cs.get_parameter()

en = time.time()

print("Time required for processing 201 waveforms using C++ CosFit: {0} s".format(en-st))