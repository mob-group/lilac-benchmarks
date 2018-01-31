import numpy as np

ref = np.array([[0., 0., 0.],   # O
                [0.0, 0.756950327264, 0.585882276618],  # H1
                [0.0, -0.756950327264, 0.585882276618], # H2
                [0.0, 0.0, 0.15]]) # M

def rotmat(p):
    theta2 = np.dot(p,p)

    if (theta2<1.0e-12):  # Then use the small-rotation formula
        RM = np.eye(3)    # Start from the identity matrix
        RM[0][1] = -p[2]  # and add first-order corrections
        RM[1][0] = p[2]
        RM[0][2] = p[1]
        RM[2][0] = -p[1]
        RM[1][2] = -p[0]
        RM[2][1] = p[0]
        return RM

    theta = np.sqrt(theta2)
    cost = np.cos(theta)
    sint = np.sin(theta)

    pn = p/theta # Normalised

    E = np.zeros((3,3))  # Skew-symmetric matrix
    E[0][1] =  -pn[2]
    E[0][2]  =  pn[1]
    E[1][2]  = -pn[0]
    E[1][0]  = -E[0][1]
    E[2][0]  = -E[0][2]
    E[2][1]  = -E[1][2]    
    Esq = np.dot(E,E)

    # Rotation matrix from Rodrigues' formula
    RM = np.eye(3)  # Rotation matrix (initially set to the identity)
    for i in xrange(3):
        for j in xrange(3):
            RM[i][j] += (1.0-cost)*Esq[i][j] + sint*E[i][j]

    return RM

def to_atomistic(rb):
    """ rb should be coords in CoM+AA form, with shape (-1,3)"""
    nmol = rb.shape[0]/2
    coords = np.zeros((4*nmol,3))

    for mol in xrange(nmol):
        CoM = rb[mol]
        p = rb[nmol+mol]

        rmat = rotmat(p)

        for i in xrange(4):
            coords[4*mol+i] = CoM + np.dot(rmat, ref[i])

    return coords
