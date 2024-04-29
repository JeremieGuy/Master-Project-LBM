import numpy as np
import sys
import matplotlib.pyplot as plt

nx = 300
ny = 100

def testMeshGrid():
    nx = np.linspace(-5,5,20)
    ny = np.linspace(-5,5,20)
    X,Y = np.meshgrid(nx,ny)

    # print(X)

    plt.plot(X,Y, 'o')
    plt.show()

def testMatrix():
    x = np.linspace(0,300)
    
    y = np.linspace(0,100,100)
    # print(y)
    # plt.plot(x,y, 'o')
    # plt.show()
    print(np.ones((30,10)))

def testFull():
    x = np.full((10,2),False)
    print(x)

def testDivide():
    x = 100
    print(x/4)
    print(int(x/4))
    print(x//4)

def testRange():
    for i in range(0,14):
        print(i)
    
def rolling():
    x =  [  [[ 0, 1,  2,  3], 
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11]],
            
            [[ 0, 1,  2,  3], 
            [ 4,  5,  6,  7],
            [ 8,  9, 10, 11]]]
    
def testSum():
    x = np.ones([10,5,3])
    b = [0,-1,1]
    print(x)
    print(sum(x,2))
    print(sum(x*b,2))

def testLinesapce():
    # x = np.arange(0,100,1)
    # print(x[1:49])
    # print(len(x[1:49]))
    x = np.arange(0,100,1)
    print(x)
    print(len(x))

def testVel(nx,ny):
    np.set_printoptions(threshold=sys.maxsize)
    lattice = np.zeros((2,nx,ny))
    lattice[0,0,:] = 0.1
    print(np.size(lattice))
    print(lattice[0,:,:])
    plt.imshow(np.transpose(lattice[0,:,:]))
    plt.show()

def testUmax():
    distance = abs(np.arange(-ny//2,ny//2,1))
    print(distance)
    print(len(distance))
    umax = 0.04
    expectedU = [umax*(1-(i/ny)**2) for i in distance]
    plt.plot(np.arange(0,ny,1),expectedU)
    plt.show()

def testcoord(x):
    distance = abs(np.arange(-ny//2,ny//2,1))
    print(distance[int(x)])

# testMatrix()
# testMeshGrid()
# testFull()
# testDivide()
# testRange()
# rolling()
# testSum()
# testLinesapce()
# testVel(30,10)
# testUmax()
testcoord("20:30")