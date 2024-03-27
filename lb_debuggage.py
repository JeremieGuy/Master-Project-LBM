import numpy as np

def distance(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

def main():
    Nx = 400
    Ny = 100
    tau = .53
    Nt = 3000

    NL = 9
    cxs = np.array([0,0,1,1,1,0,-1,-1,-1])
    cys = np.array([0,1,1,0,-1,-1,-1,0,1])
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])

    F = np.ones((Ny,Nx,NL)) + .01 * np.random.randn(Ny, Nx, NL)
    F[:, :, 3] = 2.3

    cylinder = np.full((Ny, Nx), False)

    for y in range(0,Ny):
        for x in range(0,Nx):
            if(distance(Nx//4, Ny//2, x, y)<13):
                cylinder[y][x] = True
    
    for it in range(Nt):
        print(it)

        for i, cx, cy in zip(range(NL), cxs, cys):
            F[:, :, i] = np.roll(F[:, :, i], cx, axis = 1)
            F[:, :, i] = np.roll(F[:, :, i], cy, axis = 0)
        
        bndryF = F[cylinder, :]
        bndryF = bndryF[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]]

        rho = np.sum(F, 2)
        ux = np.sum(F * cxs, 2) / rho
        uy = np.sum(F * cys, 2) / rho

        F[cylinder, :] = bndryF
        ux[cylinder] = 0
        uy[cylinder] = 0

        Feq = np.zeros(F.shape)
        for i, cx, cy, w in zip(range(NL), cxs, cys, weights):
            Feq[:, :, i] = rho * w * (
                1 + 3 * (cx*ux + cy*uy) + 9 * (cx*ux + cy*uy)**2 / 2 - 3 * (ux**2 + uy**2)/2
            )
        
        F = F + -(1/tau) * (F - Feq)


main()