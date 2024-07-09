import numpy as np
from matplotlib import pyplot as plt 
from pylab import *

nx = 300
ny = 200
# cx, cy, r = Nx//4, Ny//2, Ny//9 # Coordinates of the cylinder.

# obstacle = np.full((ny, nx), False)

# for y in range(0,ny):
#         for x in range(0,nx):
#                 if y>50 and y<150 and x<150:
#                         obstacle[y][x] = True
#                 if y>50 and y<150 and x>200 and x<250:
#                         obstacle[y][x] = True

def obstacle_fun(x, y):
#     y1 = y > 50
#     y2 = y < 150
#     x1 = x < 150
#     x2 = x > 200
#     x3 = x < 250

#     return (y1 & y2 & x1) | (y1 & y2 & x2 & x3)
        y1 = y >= 35
        y2 = y <= 65
        x1 = x >= 135
        x2 = x <= 165

        return (y1 & y2 & x1 & x2)

obstacle = fromfunction(obstacle_fun, (nx,ny))

def obstacle_fun2(x, y):
    cx = 150
    cy = 100
    r = 20
    return (x-cx)**2+(y-cy)**2<r**2

obstacle2 = fromfunction(obstacle_fun, (nx,ny))

def obstacle_fun3(x, y):
    x1 = x >= 50
    x2 = x <= 150
    y1 = y <= 150
    y2 = y >= 200
    y3 = y <= 250

    return (x1 & x2 & y1) | (x1 & x2 & y2 & y3)

obstacle3 = fromfunction(obstacle_fun3, (ny,nx))

flags = np.zeros((ny,nx))
flags[obstacle3] = 1


# figure(1)
# plt.imshow(np.transpose(obstacle3), interpolation='nearest')
plt.imshow(obstacle3)
# plt.plot(obstacle2)
plt.title("obstacle")
plt.show()

# print(obstacle2)

plt.figure
plt.imshow(flags)
plt.show()