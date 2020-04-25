import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import math

t = np.arange(0.0, 2.2, 0.001)
y1 = list(t)
x1 = list(map(lambda y: math.cos(y) + 1.0, y1))  
x2 = list(t)
y2 = list(map(lambda x: math.log10(x + 1.0) + 1.0, x2)) 

fig = plt.figure()
ax1 = fig.add_subplot(111)

line1, = ax1.plot(x1, y1, 'g', label="x1 - cos(x2) = 1")
ax1.set_xlabel('x1')
ax1.set_ylabel('x2', color='g')


# create line plot of y2(x)
line2, = ax1.plot(x2, y2, 'r', label="x2 - lg(x1 + 1) = 1")

# set title, plot limits, etc
plt.title('System illustrate:')
plt.xlim(0.0, 2.0)
plt.ylim(0.0, 2.0)

# add a legend, and position it on the upper right
plt.legend((line1, line2), ("x1 - cos(x2) = 1", "x2 - lg(x1 + 1) = 1"))

plt.show()